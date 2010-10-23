/**
 * copyright:
 * Research Group on Numerical Methods for PDEs
 * University of Göttingen
 *
 * lube@math.uni-goettingen.de
 */

#ifndef TIMEBLOCK_H
#define TIMEBLOCK_H

template<typename STREAM>
class TimeBlock
{
  public:
    TimeBlock(STREAM & stream_, const char* blockname, bool singleline=true)
		    : m_blockname(blockname),
		      m_singleline(singleline),
		      m_running(true),
		      m_timer(MPI_COMM_WORLD, true),
		      stream(stream_)
      {
	MPI_Barrier(MPI_COMM_WORLD);
	stream << blockname << " ... ";
	if (singleline)
	  stream << std::flush;
	else
	  stream << std::endl;
	m_timer.start();
      }
    
    ~TimeBlock()
      {
	close();
      }

    void close()
      {
	if (!m_running)
	  return;
	m_running = false;
	
	m_timer.stop();
	
	if (!m_singleline)
	  stream << m_blockname << " took ";

	stream << std::fixed << std::setprecision(4);
	
	m_timer.print_data(stream);
	
	MPI_Barrier(MPI_COMM_WORLD);
      }


  private:
    const char* m_blockname;
    bool m_singleline;
    bool m_running;
    dealii::Timer m_timer;
    STREAM & stream;
};

#endif

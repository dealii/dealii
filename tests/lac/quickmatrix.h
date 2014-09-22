
/**
 * Hard-coded Laplacian matrix.
 * Just a quick matrix to investigate processor performance.
 * It implements a finite difference scheme on a grid of grid size 1
 * with #nx# times #ny# grid points.
 * The diagonal is scaled to 1, resulting in an effective mesh width of 1/2.
 */

template<typename number>
class QuickMatrix
{
public:
  /**
   * Constructor initializing the grid size.
   */
  QuickMatrix(unsigned int nx, unsigned int ny);

  /**
   * Matrix-vector-product.
   */
  template <typename number2>
  void vmult(Vector<number2> &, const Vector<number2> &) const;
protected:
  const unsigned int nx;
  const unsigned int ny;
};


template<typename number>
QuickMatrix<number>::QuickMatrix(unsigned int nx, unsigned int ny)
  :
  nx(nx), ny(ny)
{}

template<typename number>
template <typename number2>
void
QuickMatrix<number>::vmult(Vector<number2> &d,
                           const Vector<number2> &s) const
{
  const unsigned int step = nx-1;
  const unsigned int right = step-1;
  const unsigned int top = ny-1;

  // Bottom row

  d(0) = s(0) - .25 * ( s(1) + s(step) );

  for (unsigned int x=1; x<right; ++x)
    d(x) = s(x) - .25 * ( s(x-1) + s(x+1) + s(x+step) );

  d(right) = s(right) - .25 * ( s(right-1) + s(right+step) );

  // Middle rows

  unsigned int start = 0;
  for (unsigned int y=1; y<top; y++)
    {
      start += step;
      d(start) = s(start)
                 - .25 * ( s(start-step) + s(start+1) + s(start+step) );

      for (unsigned int x=1; x<right; ++x)
        {
          const unsigned int xy = start+x;
          d(xy) = s(xy)
                  - .25 * ( s(xy-step) + s(xy-1) + s(xy+1) + s(xy+step) );
        }
      d(start+right) = s(start+right)
                       - .25 * ( s(start+right-1) + s(start+right+step) );
    }

  // Top row

  start += step;
  d(start) = s(start)
             - .25 * ( s(start-step) + s(start+1) );

  for (unsigned int x=1; x<right; ++x)
    {
      const unsigned int xy = start+x;
      d(xy) = s(xy)
              - .25 * ( s(xy-step) + s(xy-1) + s(xy+1) );
    }
  d(start+right) = s(start+right)
                   - .25 * ( s(start+right-step) + s(start+right-1) );
}

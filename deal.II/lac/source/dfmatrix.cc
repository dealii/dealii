// $Id$

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <lac/dfmatrix.h>





dFMatrix::dFMatrix (const dFMatrix &m)
{
  init (m.dim_image, m.dim_range);
  for (int i=0; i!=dim_image*dim_range; ++i)
    val[i] = m.val[i];
};

void dFMatrix::init(int mm, int nn)
{
  val_size = nn*mm;
  val = new double[val_size];
  dim_range = nn;
  dim_image = mm;
  clear ();
}

dFMatrix::~dFMatrix()
{
  delete[] val;
}

void dFMatrix::reinit(int mm, int nn)
{
  if (val_size<nn*mm)
    {
      delete[] val;
      init(mm, nn);
    }
  else
    {
      dim_range = nn;
      dim_image = mm;
//      memset(val, 0, sizeof(double)*nn*mm);
      for (int i = nn*mm-1; i>=0 ; i--) val[i] = 0;
    }
}

void dFMatrix::vmult(dVector& dst, const dVector& src,const int adding) const
{
  Assert(dst.n() == m(), ExcDimensionMismatch(dst.n(), m()));
  Assert(src.n() == n(), ExcDimensionMismatch(src.n(), n()));

  double s;
  if ((n()==3) && (m()==3))
  {
    double s0,s1,s2;
    s   = src(0);
    s0  = s*val[0]; s1  = s*val[3]; s2  = s*val[6]; 
    s   = src(1);
    s0 += s*val[1]; s1 += s*val[4]; s2 += s*val[7];
    s   = src(2);
    s0 += s*val[2]; s1 += s*val[5]; s2 += s*val[8];

    if (!adding)
    {
      dst(0) = s0;
      dst(1) = s1;
      dst(2) = s2;
    }
    else
    {
      dst(0) += s0;
      dst(1) += s1;
      dst(2) += s2;
    }
  }
  else if ((n()==4) && (m()==4))
  {
    double s0,s1,s2,s3;
    s = src(0);
    s0 = s*val[0]; s1 = s*val[4]; s2 = s*val[8]; s3 = s*val[12];
    s = src(1);
    s0 += s*val[1]; s1 += s*val[5]; s2 += s*val[9]; s3 += s*val[13];
    s = src(2);
    s0 += s*val[2]; s1 += s*val[6]; s2 += s*val[10]; s3 += s*val[14];
    s = src(3);
    s0 += s*val[3]; s1 += s*val[7]; s2 += s*val[11]; s3 += s*val[15];
    
    if (!adding)
    {
      dst(0) = s0;
      dst(1) = s1;
      dst(2) = s2;
      dst(3) = s3;
    }
    else
    {
      dst(0) += s0;
      dst(1) += s1;
      dst(2) += s2;
      dst(3) += s3;
    }
  }
  else if ((n()==8) && (m()==8))
  {
    double s0,s1,s2,s3,s4,s5,s6,s7;
    s = src(0);
    s0 = s*val[0]; s1 = s*val[8]; s2 = s*val[16]; s3 = s*val[24];
    s4 = s*val[32]; s5 = s*val[40]; s6 = s*val[48]; s7 = s*val[56];
    s = src(1);
    s0 += s*val[1]; s1 += s*val[9]; s2 += s*val[17]; s3 += s*val[25];
    s4 += s*val[33]; s5 += s*val[41]; s6 += s*val[49]; s7 += s*val[57];
    s = src(2);
    s0 += s*val[2]; s1 += s*val[10]; s2 += s*val[18]; s3 += s*val[26];
    s4 += s*val[34]; s5 += s*val[42]; s6 += s*val[50]; s7 += s*val[58];
    s = src(3);
    s0 += s*val[3]; s1 += s*val[11]; s2 += s*val[19]; s3 += s*val[27];
    s4 += s*val[35]; s5 += s*val[43]; s6 += s*val[51]; s7 += s*val[59];
    s = src(4);
    s0 += s*val[4]; s1 += s*val[12]; s2 += s*val[20]; s3 += s*val[28];
    s4 += s*val[36]; s5 += s*val[44]; s6 += s*val[52]; s7 += s*val[60];
    s = src(5);
    s0 += s*val[5]; s1 += s*val[13]; s2 += s*val[21]; s3 += s*val[29];
    s4 += s*val[37]; s5 += s*val[45]; s6 += s*val[53]; s7 += s*val[61];
    s = src(6);
    s0 += s*val[6]; s1 += s*val[14]; s2 += s*val[22]; s3 += s*val[30];
    s4 += s*val[38]; s5 += s*val[46]; s6 += s*val[54]; s7 += s*val[62];
    s = src(7);
    s0 += s*val[7]; s1 += s*val[15]; s2 += s*val[23]; s3 += s*val[31];
    s4 += s*val[39]; s5 += s*val[47]; s6 += s*val[55]; s7 += s*val[63];
    
    if (!adding)
    {
      dst(0) = s0;
      dst(1) = s1;
      dst(2) = s2;
      dst(3) = s3;
      dst(4) = s4;
      dst(5) = s5;
      dst(6) = s6;
      dst(7) = s7;
    }
    else
    {
      dst(0) += s0;
      dst(1) += s1;
      dst(2) += s2;
      dst(3) += s3;
      dst(4) += s4;
      dst(5) += s5;
      dst(6) += s6;
      dst(7) += s7;
    }
  }
  else
  {    
    double* e = val;
    for (int i=0;i<m();i++)
    {
      s = 0.;
      for (int j=0;j<n();j++) s += src(j) * *(e++);
      if (!adding) dst(i) = s;
      else dst(i) += s;
    }
  }
}

void dFMatrix::gsmult(dVector& dst, const dVector& src,const iVector& gl) const
{
  Assert(n() == m(), ExcNotQuadratic());
  Assert(dst.n() == n(), ExcDimensionMismatch(dst.n(), n()));
  Assert(src.n() == n(), ExcDimensionMismatch(src.n(), n()));
  Assert(gl.n() == n(), ExcDimensionMismatch(gl.n(), n()));

  double s;
  if ((n()==3) && (m()==3))
  {
    double s0=0.,s1=0.,s2=0.;
    s   = src(0);
    if(gl(1)<gl(0)) s1  = s*val[3]; if(gl(2)<gl(0))  s2  = s*val[6]; 
    s   = src(1);
    if(gl(0)<gl(1)) s0 += s*val[1]; if(gl(2)<gl(1))  s2 += s*val[7];
    s   = src(2);
    if(gl(0)<gl(2)) s0 += s*val[2]; if(gl(1)<gl(2))  s1 += s*val[5]; 

    dst(0) += s0;
    dst(1) += s1;
    dst(2) += s2;
  }
  else if ((n()==4) && (m()==4))
  {
    double s0=0.,s1=0.,s2=0.,s3=0.;
    s = src(0);
    if(gl(1)<gl(0)) s1 = s*val[4];  if(gl(2)<gl(0)) s2 = s*val[8]; if(gl(3)<gl(0)) s3 = s*val[12];
    s = src(1);
    if(gl(0)<gl(1)) s0 += s*val[1]; if(gl(2)<gl(1)) s2 += s*val[9]; if(gl(3)<gl(1)) s3 += s*val[13];
    s = src(2);
    if(gl(0)<gl(2)) s0 += s*val[2]; if(gl(1)<gl(2)) s1 += s*val[6]; if(gl(3)<gl(2)) s3 += s*val[14];
    s = src(3);
    if(gl(0)<gl(3)) s0 += s*val[3]; if(gl(1)<gl(3)) s1 += s*val[7]; if(gl(2)<gl(3)) s2 += s*val[11];

    dst(0) += s0;
    dst(1) += s1;
    dst(2) += s2;
    dst(3) += s3;
  }
  else if ((n()==8) && (m()==8))
  {
    double s0=0.,s1=0.,s2=0.,s3=0.,s4=0.,s5=0.,s6=0.,s7=0.;
    s = src(0);
    if(gl(1)<gl(0)) s1 = s*val[8]; 
    if(gl(2)<gl(0)) s2 = s*val[16]; 
    if(gl(3)<gl(0)) s3 = s*val[24];
    if(gl(4)<gl(0)) s4 = s*val[32]; 
    if(gl(5)<gl(0)) s5 = s*val[40]; 
    if(gl(6)<gl(0)) s6 = s*val[48]; 
    if(gl(7)<gl(0)) s7 = s*val[56];
    s = src(1);
    if(gl(0)<gl(1)) s0 += s*val[1]; 
    if(gl(2)<gl(1)) s2 += s*val[17]; 
    if(gl(3)<gl(1)) s3 += s*val[25];
    if(gl(4)<gl(1)) s4 += s*val[33]; 
    if(gl(5)<gl(1)) s5 += s*val[41]; 
    if(gl(6)<gl(1)) s6 += s*val[49]; 
    if(gl(7)<gl(1)) s7 += s*val[57];
    s = src(2);
    if(gl(0)<gl(2)) s0 += s*val[2]; 
    if(gl(1)<gl(2)) s1 += s*val[10]; 
    if(gl(3)<gl(2)) s3 += s*val[26];
    if(gl(4)<gl(2)) s4 += s*val[34]; 
    if(gl(5)<gl(2)) s5 += s*val[42]; 
    if(gl(6)<gl(2)) s6 += s*val[50]; 
    if(gl(7)<gl(2)) s7 += s*val[58];
    s = src(3);
    if(gl(0)<gl(3)) s0 += s*val[3]; 
    if(gl(1)<gl(3)) s1 += s*val[11]; 
    if(gl(2)<gl(3)) s2 += s*val[19]; 
    if(gl(4)<gl(3)) s4 += s*val[35]; 
    if(gl(5)<gl(3)) s5 += s*val[43]; 
    if(gl(6)<gl(3)) s6 += s*val[51]; 
    if(gl(7)<gl(3)) s7 += s*val[59];
    s = src(4);
    if(gl(0)<gl(4)) s0 += s*val[4]; 
    if(gl(1)<gl(4)) s1 += s*val[12]; 
    if(gl(2)<gl(4)) s2 += s*val[20]; 
    if(gl(3)<gl(4)) s3 += s*val[28];
    if(gl(5)<gl(4)) s5 += s*val[44]; 
    if(gl(6)<gl(4)) s6 += s*val[52]; 
    if(gl(7)<gl(4)) s7 += s*val[60];
    s = src(5);
    if(gl(0)<gl(5)) s0 += s*val[5]; 
    if(gl(1)<gl(5)) s1 += s*val[13]; 
    if(gl(2)<gl(5)) s2 += s*val[21]; 
    if(gl(3)<gl(5)) s3 += s*val[29];
    if(gl(4)<gl(5)) s4 += s*val[37]; 
    if(gl(6)<gl(5)) s6 += s*val[53]; 
    if(gl(7)<gl(5)) s7 += s*val[61];
    s = src(6);
    if(gl(0)<gl(6)) s0 += s*val[6]; 
    if(gl(1)<gl(6)) s1 += s*val[14]; 
    if(gl(2)<gl(6)) s2 += s*val[22]; 
    if(gl(3)<gl(6)) s3 += s*val[30];
    if(gl(4)<gl(6)) s4 += s*val[38]; 
    if(gl(5)<gl(6)) s5 += s*val[46]; 
    if(gl(7)<gl(6)) s7 += s*val[62];
    s = src(7);
    if(gl(0)<gl(7)) s0 += s*val[7]; 
    if(gl(1)<gl(7)) s1 += s*val[15]; 
    if(gl(2)<gl(7)) s2 += s*val[23]; 
    if(gl(3)<gl(7)) s3 += s*val[31];
    if(gl(4)<gl(7)) s4 += s*val[39]; 
    if(gl(5)<gl(7)) s5 += s*val[47]; 
    if(gl(6)<gl(7)) s6 += s*val[55]; 
    
    dst(0) += s0;
    dst(1) += s1;
    dst(2) += s2;
    dst(3) += s3;
    dst(4) += s4;
    dst(5) += s5;
    dst(6) += s6;
    dst(7) += s7;
  }
  else
  {    
    double* e = val;
    for (int i=0;i<m();i++)
    {
      s = 0.;
      for (int j=0;j<n();j++) if(gl(i)<gl(j)) s += src(j) * *(e++);
      dst(i) += s;
    }
  }
}

void dFMatrix::Tvmult(dVector& dst, const dVector& src, const int adding) const
{
  Assert(dst.n() == n(), ExcDimensionMismatch(dst.n(), n()));
  Assert(src.n() == m(), ExcDimensionMismatch(src.n(), m()));

  int i,j;
  double s;
  for (i=0;i<m();i++)
  {
    s = 0.;
    for (j=0;j<n();j++) s += src(j) * el(j,i);
    if(!adding) dst(i) = s;
    else dst(i) += s;
  }
}

double dFMatrix::residual(dVector& dst, const dVector& src,
			  const dVector& right) const
{
  Assert(dst.n() == m(), ExcDimensionMismatch(dst.n(), m()));
  Assert(src.n() == n(), ExcDimensionMismatch(src.n(), n()));
  Assert(right.n() == m(), ExcDimensionMismatch(right.n(), m()));

  int i,j;
  double s, res = 0.;
  for (i=0;i<n();i++)
    {
      s = right(i);
      for (j=0;j<m();j++) s -= src(j) * el(i,j);
      dst(i) = s;
      res += s*s;
    }
  return sqrt(res);
}

void dFMatrix::forward(dVector& dst, const dVector& src) const
{
  Assert(n() == m(), ExcNotQuadratic());
  Assert(dst.n() == n(), ExcDimensionMismatch(dst.n(), n()));
  Assert(src.n() == n(), ExcDimensionMismatch(src.n(), n()));

  int i,j;
  int nu = MIN(m(),n());
  double s;
  for (i=0;i<nu;i++)
    {
      s = src(i);
      for (j=0;j<i;j++) s -= dst(j) * el(i,j);
      dst(i) = s/el(i,i);
    }
}

void dFMatrix::backward(dVector& dst, const dVector& src) const
{
  Assert(n() == m(), ExcNotQuadratic());
  Assert(dst.n() == n(), ExcDimensionMismatch(dst.n(), n()));
  Assert(src.n() == n(), ExcDimensionMismatch(src.n(), n()));

  int i,j;
  int nu = MIN(m(),n());
  double s;
  for (i=nu-1;i>=0;i--)
    {
      s = src(i);
      for (j=i+1;j<nu;j++) s -= dst(j) * el(i,j);
      dst(i) = s/el(i,i);
    }
}

dFMatrix& dFMatrix::operator = (const dFMatrix& M)
{
  reinit(M);
  int nn = n()*m();
  for (int i=0;i<nn;i++) val[i] = M.val[i];
  return *this;
}

void dFMatrix::fill(const dFMatrix& src, int i, int j)
{
  Assert (n() >= src.n() + j, ExcInvalidDestination(n(), src.n(), j));
  Assert (m() >= src.m() + i, ExcInvalidDestination(m(), src.m(), i));

  for (int ii=0; ii<src.m() ; ii++)
    for (int jj=0; jj<src.n() ; jj++)
      el(ii+i,jj+j) = src.el(ii,jj);
}

void dFMatrix::add_row(int i, double s, int j)
{
  int k;
  for (k=0;k<m();k++) el(i,k) += s*el(j,k);
}

void dFMatrix::add_row(int i, double s, int j, double t, int k)
{
  int l;
  for (l=0;l<m();l++) el(i,l) += s*el(j,l) + t*el(k,l);
}

void dFMatrix::add_col(int i, double s, int j)
{
  int k;
  for (k=0;k<n();k++) el(k,i) += s*el(k,j);
}

void dFMatrix::add_col(int i, double s, int j, double t, int k)
{
  int l;
  for (l=0;l<n();l++) el(l,i) += s*el(l,j) + t*el(l,k);
}

void dFMatrix::swap_row(int i, int j)
{
  int k;
  double s;
  for (k=0;k<m();k++)
  {
    s = el(i,k); el(i,k) = el(j,k); el(j,k) = s;
  }
}

void dFMatrix::swap_col(int i, int j)
{
  int k;
  double s;
  for (k=0;k<n();k++)
  {
    s = el(k,i); el(k,i) = el(k,j); el(k,j) = s;
  }
}

void dFMatrix::diagadd(const double& src)
{
  Assert (m() == n(), ExcDimensionMismatch(m(),n()));
  for (int i=0;i<n();i++)
    el(i,i) += src;
}

void dFMatrix::mmult(dFMatrix& dst, const dFMatrix& src) const
{
  Assert (n() == src.m(), ExcDimensionMismatch(n(), src.m()));
  int i,j,k;
  double s = 1.;
  dst.reinit(m(), src.n());

  for (i=0;i<m();i++)
    for (j=0;j<src.n();j++)
      {
	s = 0.;
	for (k=0;k<n();k++) s+= el(i,k) * src.el(k,j);
	dst.el(i,j) = s;
      }
}

/*void dFMatrix::mmult(dFMatrix& dst, const dFMatrix& src) const
{
  Assert (m() == src.n(), ExcDimensionMismatch(m(), src.n()));

  int i,j,k;
  double s = 1.;

  dst.reinit(n(), src.m());

  for (i=0;i<n();i++)
    for (j=0;j<src.m();j++)
      {
	s = 0.;
	for (k=0;k<m();k++) s+= el(i,k) * src.el(k,j);
	dst.el(i,j) = s;
      }
}*/

void dFMatrix::Tmmult(dFMatrix& dst, const dFMatrix& src) const
{
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));

  int i,j,k;
  double s = 1.;
  dst.reinit(m(), src.m());

  for (i=0;i<m();i++)
    for (j=0;j<src.m();j++)
      {
	s = 0.;
	for (k=0;k<n();k++) s+= el(k,i) * src.el(k,j);
	dst.el(i,j) = s;
      }
}

/*void dFMatrix::Tmmult(dFMatrix& dst, const dFMatrix& src) const
{
  Assert (m() == src.n(), ExcDimensionMismatch(m(), src.n()));

  int i,j,k;
  double s = 1.;
  
  dst.reinit(n(), src.m());

  for (i=0;i<n();i++)
    for (j=0;j<src.m();j++)
      {
	s = 0.;
	for (k=0;k<m();k++) s+= el(k,i) * src.el(k,j);
	dst.el(i,j) = s;
      }
}*/

void dFMatrix::print(FILE* f, const char* format) const
{
  if (!format) format = " %5.2f";
  int i,j;
  for (i=0;i<m();i++)
    {
      for (j=0;j<n();j++) fprintf(f, format, el(i,j));
      fprintf(f, "\n");
    }
}

void dFMatrix::add(double s,const dFMatrix& src)
{
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));
  if ((n()==3) && (m()==3))
  {
    val[0] += s * src.el(0);
    val[1] += s * src.el(1);
    val[2] += s * src.el(2);
    val[3] += s * src.el(3);
    val[4] += s * src.el(4);
    val[5] += s * src.el(5);
    val[6] += s * src.el(6);
    val[7] += s * src.el(7);
    val[8] += s * src.el(8);
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * src.el(0);
    val[1] += s * src.el(1);
    val[2] += s * src.el(2);
    val[3] += s * src.el(3);
    val[4] += s * src.el(4);
    val[5] += s * src.el(5);
    val[6] += s * src.el(6);
    val[7] += s * src.el(7);
    val[8] += s * src.el(8);
    val[9] += s * src.el(9);
    val[10] += s * src.el(10);
    val[11] += s * src.el(11);
    val[12] += s * src.el(12);
    val[13] += s * src.el(13);
    val[14] += s * src.el(14);
    val[15] += s * src.el(15);
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * src.el(0);
    val[1] += s * src.el(1);
    val[2] += s * src.el(2);
    val[3] += s * src.el(3);
    val[4] += s * src.el(4);
    val[5] += s * src.el(5);
    val[6] += s * src.el(6);
    val[7] += s * src.el(7);
    val[8] += s * src.el(8);
    val[9] += s * src.el(9);
    val[10] += s * src.el(10);
    val[11] += s * src.el(11);
    val[12] += s * src.el(12);
    val[13] += s * src.el(13);
    val[14] += s * src.el(14);
    val[15] += s * src.el(15);
    val[16] += s * src.el(16);
    val[17] += s * src.el(17);
    val[18] += s * src.el(18);
    val[19] += s * src.el(19);

    val[20] += s * src.el(20);
    val[21] += s * src.el(21);
    val[22] += s * src.el(22);
    val[23] += s * src.el(23);
    val[24] += s * src.el(24);
    val[25] += s * src.el(25);
    val[26] += s * src.el(26);
    val[27] += s * src.el(27);
    val[28] += s * src.el(28);
    val[29] += s * src.el(29);

    val[30] += s * src.el(30);
    val[31] += s * src.el(31);
    val[32] += s * src.el(32);
    val[33] += s * src.el(33);
    val[34] += s * src.el(34);
    val[35] += s * src.el(35);
    val[36] += s * src.el(36);
    val[37] += s * src.el(37);
    val[38] += s * src.el(38);
    val[39] += s * src.el(39);

    val[40] += s * src.el(40);
    val[41] += s * src.el(41);
    val[42] += s * src.el(42);
    val[43] += s * src.el(43);
    val[44] += s * src.el(44);
    val[45] += s * src.el(45);
    val[46] += s * src.el(46);
    val[47] += s * src.el(47);
    val[48] += s * src.el(48);
    val[49] += s * src.el(49);

    val[50] += s * src.el(50);
    val[51] += s * src.el(51);
    val[52] += s * src.el(52);
    val[53] += s * src.el(53);
    val[54] += s * src.el(54);
    val[55] += s * src.el(55);
    val[56] += s * src.el(56);
    val[57] += s * src.el(57);
    val[58] += s * src.el(58);
    val[59] += s * src.el(59);

    val[60] += s * src.el(60);
    val[61] += s * src.el(61);
    val[62] += s * src.el(62);
    val[63] += s * src.el(63);
  }
  else
  {
    for (int i = n()*m()-1; i>=0 ; i--) val[i] += s * src.el(i);
  }
}



void dFMatrix::add_diag(double s,const dFMatrix& src)
{
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));

  if ((n()==3) && (m()==3))
  {
    val[0] += s * src.el(0);
    val[0] += s * src.el(1);
    val[0] += s * src.el(2);
    val[3] += s * src.el(3);
    val[3] += s * src.el(4);
    val[3] += s * src.el(5);
    val[6] += s * src.el(6);
    val[6] += s * src.el(7);
    val[6] += s * src.el(8);
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * src.el(0);
    val[0] += s * src.el(1);
    val[0] += s * src.el(2);
    val[0] += s * src.el(3);
    val[4] += s * src.el(4);
    val[4] += s * src.el(5);
    val[4] += s * src.el(6);
    val[4] += s * src.el(7);
    val[8] += s * src.el(8);
    val[8] += s * src.el(9);
    val[8] += s * src.el(10);
    val[8] += s * src.el(11);
    val[12] += s * src.el(12);
    val[12] += s * src.el(13);
    val[12] += s * src.el(14);
    val[12] += s * src.el(15);
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * src.el(0);
    val[0] += s * src.el(1);
    val[0] += s * src.el(2);
    val[0] += s * src.el(3);
    val[0] += s * src.el(4);
    val[0] += s * src.el(5);
    val[0] += s * src.el(6);
    val[0] += s * src.el(7);
    val[8] += s * src.el(8);
    val[8] += s * src.el(9);
    val[8] += s * src.el(10);
    val[8] += s * src.el(11);
    val[8] += s * src.el(12);
    val[8] += s * src.el(13);
    val[8] += s * src.el(14);
    val[8] += s * src.el(15);
    val[16] += s * src.el(16);
    val[16] += s * src.el(17);
    val[16] += s * src.el(18);
    val[16] += s * src.el(19);

    val[16] += s * src.el(20);
    val[16] += s * src.el(21);
    val[16] += s * src.el(22);
    val[16] += s * src.el(23);
    val[24] += s * src.el(24);
    val[24] += s * src.el(25);
    val[24] += s * src.el(26);
    val[24] += s * src.el(27);
    val[24] += s * src.el(28);
    val[24] += s * src.el(29);

    val[24] += s * src.el(30);
    val[24] += s * src.el(31);
    val[32] += s * src.el(32);
    val[32] += s * src.el(33);
    val[32] += s * src.el(34);
    val[32] += s * src.el(35);
    val[32] += s * src.el(36);
    val[32] += s * src.el(37);
    val[32] += s * src.el(38);
    val[32] += s * src.el(39);

    val[40] += s * src.el(40);
    val[40] += s * src.el(41);
    val[40] += s * src.el(42);
    val[40] += s * src.el(43);
    val[40] += s * src.el(44);
    val[40] += s * src.el(45);
    val[40] += s * src.el(46);
    val[40] += s * src.el(47);
    val[48] += s * src.el(48);
    val[48] += s * src.el(49);

    val[48] += s * src.el(50);
    val[48] += s * src.el(51);
    val[48] += s * src.el(52);
    val[48] += s * src.el(53);
    val[48] += s * src.el(54);
    val[48] += s * src.el(55);
    val[56] += s * src.el(56);
    val[56] += s * src.el(57);
    val[56] += s * src.el(58);
    val[56] += s * src.el(59);

    val[56] += s * src.el(60);
    val[56] += s * src.el(61);
    val[56] += s * src.el(62);
    val[56] += s * src.el(63);
  }
  else
  {
    for (int i = n()*m()-1; i>=0 ; i--) val[i] += s * src.el(i);
  }
}

void dFMatrix::Tadd(double s,const dFMatrix& src)
{
  Assert (m() == n(),     ExcNotQuadratic());
  Assert (m() == src.m(), ExcDimensionMismatch(m(), src.m()));
  Assert (n() == src.n(), ExcDimensionMismatch(n(), src.n()));

  if ((n()==3) && (m()==3))
  {
    val[0] += s * src.el(0);
    val[1] += s * src.el(3);
    val[2] += s * src.el(6);

    val[3] += s * src.el(1);
    val[4] += s * src.el(4);
    val[5] += s * src.el(7);

    val[6] += s * src.el(2);
    val[7] += s * src.el(5);
    val[8] += s * src.el(8);
  }
  else if ((n()==4) && (m()==4))
  {
    val[0] += s * src.el(0);
    val[1] += s * src.el(4);
    val[2] += s * src.el(8);
    val[3] += s * src.el(12);

    val[4] += s * src.el(1);
    val[5] += s * src.el(5);
    val[6] += s * src.el(9);
    val[7] += s * src.el(13);

    val[8] += s * src.el(2);
    val[9] += s * src.el(6);
    val[10] += s * src.el(10);
    val[11] += s * src.el(14);

    val[12] += s * src.el(3);
    val[13] += s * src.el(7);
    val[14] += s * src.el(11);
    val[15] += s * src.el(15);
  }
  else if ((n()==8) && (m()==8))
  {
    val[0] += s * src.el(0);
    val[1] += s * src.el(8);
    val[2] += s * src.el(16);
    val[3] += s * src.el(24);
    val[4] += s * src.el(32);
    val[5] += s * src.el(40);
    val[6] += s * src.el(48);
    val[7] += s * src.el(56);

    val[8] += s * src.el(1);
    val[9] += s * src.el(9);
    val[10] += s * src.el(17);
    val[11] += s * src.el(25);
    val[12] += s * src.el(33);
    val[13] += s * src.el(41);
    val[14] += s * src.el(49);
    val[15] += s * src.el(57);

    val[16] += s * src.el(2);
    val[17] += s * src.el(10);
    val[18] += s * src.el(18);
    val[19] += s * src.el(26);
    val[20] += s * src.el(34);
    val[21] += s * src.el(42);
    val[22] += s * src.el(50);
    val[23] += s * src.el(58);

    val[24] += s * src.el(3);
    val[25] += s * src.el(11);
    val[26] += s * src.el(19);
    val[27] += s * src.el(27);
    val[28] += s * src.el(35);
    val[29] += s * src.el(43);
    val[30] += s * src.el(51);
    val[31] += s * src.el(59);

    val[32] += s * src.el(4);
    val[33] += s * src.el(12);
    val[34] += s * src.el(20);
    val[35] += s * src.el(28);
    val[36] += s * src.el(36);
    val[37] += s * src.el(44);
    val[38] += s * src.el(52);
    val[39] += s * src.el(60);

    val[40] += s * src.el(5);
    val[41] += s * src.el(13);
    val[42] += s * src.el(21);
    val[43] += s * src.el(29);
    val[44] += s * src.el(37);
    val[45] += s * src.el(45);
    val[46] += s * src.el(53);
    val[47] += s * src.el(61);

    val[48] += s * src.el(6);
    val[49] += s * src.el(14);
    val[50] += s * src.el(22);
    val[51] += s * src.el(30);
    val[52] += s * src.el(38);
    val[53] += s * src.el(46);
    val[54] += s * src.el(54);
    val[55] += s * src.el(62);

    val[56] += s * src.el(7);
    val[57] += s * src.el(15);
    val[58] += s * src.el(23);
    val[59] += s * src.el(31);
    val[60] += s * src.el(39);
    val[61] += s * src.el(47);
    val[62] += s * src.el(55);
    val[63] += s * src.el(63);
  }
  else
    Assert (false, ExcInternalError());
}


bool
dFMatrix::operator == (const dFMatrix &m) const
{
  bool q = (dim_range==m.dim_range) && (dim_image==m.dim_image);
  if (!q) return false;

  for (int i=0; i<dim_image; ++i)
    for (int j=0; j<dim_range; ++j)
      if (el(i,j) != m.el(i,j)) return false;
  return true;
};


double dFMatrix::determinant () const {
  Assert (dim_range == dim_image,
	  ExcDimensionMismatch(dim_range, dim_image));
  Assert ((dim_range>=1) && (dim_range<=3), ExcNotImplemented(dim_range));
  
  switch (dim_range) 
    {
      case 1:
            return el(0,0);
      case 2:
            return el(0,0)*el(1,1) - el(1,0)*el(0,1);
      case 3:
            return  (el(0,0)*el(1,1)*el(2,2)
                     -el(0,0)*el(1,2)*el(2,1)
                     -el(1,0)*el(0,1)*el(2,2)
                     +el(1,0)*el(0,2)*el(2,1)
                     +el(2,0)*el(0,1)*el(1,2)
                     -el(2,0)*el(0,2)*el(1,1));
      default:
            return 0;
    };
};



void dFMatrix::clear () {
  for (int i=0; i<val_size; ++i)
    val[i] = 0.;
};



void dFMatrix::invert (const dFMatrix &M) {
  Assert (dim_range == dim_image, ExcNotQuadratic());
  Assert ((dim_range>=1) && (dim_range<=3), ExcNotImplemented(dim_range));
  Assert (dim_range == M.dim_range,
          ExcDimensionMismatch(dim_range,M.dim_range));
  Assert (dim_image == M.dim_image,
	  ExcDimensionMismatch(dim_image,M.dim_image));
  
  switch (dim_range) 
    {
      case 1:
	    val[0] = 1.0/M.val[0];
	    return;
      case 2:
					     // this is Maple output,
					     // thus a bit unstructured
      {
	    const double t4 = 1.0/(M.el(0,0)*M.el(1,1)-M.el(0,1)*M.el(1,0));
	    el(0,0) = M.el(1,1)*t4;
	    el(0,1) = -M.el(0,1)*t4;
	    el(1,0) = -M.el(1,0)*t4;
	    el(1,1) = M.el(0,0)*t4;
	    return;
      };
      
      case 3:
      {
	    const double t4 = M.el(0,0)*M.el(1,1),
			 t6 = M.el(0,0)*M.el(1,2),
			 t8 = M.el(0,1)*M.el(1,0),
			t00 = M.el(0,2)*M.el(1,0),
			t01 = M.el(0,1)*M.el(2,0),
			t04 = M.el(0,2)*M.el(2,0),
			t07 = 1.0/(t4*M.el(2,2)-t6*M.el(2,1)-t8*M.el(2,2)+
				   t00*M.el(2,1)+t01*M.el(1,2)-t04*M.el(1,1));
	    el(0,0) = (M.el(1,1)*M.el(2,2)-M.el(1,2)*M.el(2,1))*t07;
	    el(0,1) = -(M.el(0,1)*M.el(2,2)-M.el(0,2)*M.el(2,1))*t07;
	    el(0,2) = -(-M.el(0,1)*M.el(1,2)+M.el(0,2)*M.el(1,1))*t07;
	    el(1,0) = -(M.el(1,0)*M.el(2,2)-M.el(1,2)*M.el(2,0))*t07;
	    el(1,1) = (M.el(0,0)*M.el(2,2)-t04)*t07;
	    el(1,2) = -(t6-t00)*t07;
	    el(2,0) = -(-M.el(1,0)*M.el(2,1)+M.el(1,1)*M.el(2,0))*t07;
	    el(2,1) = -(M.el(0,0)*M.el(2,1)-t01)*t07;
	    el(2,2) = (t4-t8)*t07;
	    return;
      };
    };    
};
  

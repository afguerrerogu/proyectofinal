#include <iostream>
#include <cmath>
#include <vector>

double yi =2.5; // {x0, y0, z0, px , py, pz} 
double x0 =0.0;
double z0 =0.0; 
double px0 =0.0; 
double py0 =1.68884; 
double pz0 =0.2; 

const double A= 1.0;
const double C= 1.0;
const double w= 0.25;
const double a= 1.25;
const double b= 1.0;
const double c= 0.75;
//const double x[0]=2.5; 
//const double y[0]=0.0;
//const double z[0]=0.0;
//const double px[0]=0.0;
//const double py[0]=1.68884;
//const double pz[0]=0.2;

double f(double t,double  x,double y,double z,double px, double py,double pz , int id);
void rk4(double ta, double tb, double h, double x, double y, double z, double px, double py,double pz);

int main(void){

    rk4(0, 10000 , 0.1, x0,yi,z0,px0,py0,pz0);
    return 0;
}

double f(double t, double x,double y,double z,double px,double py,double pz , int id){

    if (0==id){
        return px+w*y;
    }
    if (1==id){
        return py-w*x;
    }
    if (2==id){
        return pz; 
    }
    if (3==id){
        return w*py-(2*A*x)/(pow(a,2)*(C+pow(x/a ,2)+pow(y/b ,2)+pow(z/c ,2)));
    }
    if (4==id){
        return -w*px-(2*A*y)/(pow(b,2)*(C+pow(x/a ,2)+pow(y/b ,2)+pow(z/c ,2)));
    }
    if (5==id){
        return -(2*A*z)/(pow(c,2)*(C+pow(x/a ,2)+pow(y/b ,2)+pow(z/c ,2)));
    }
    else
    {
        exit(1);
    }
}

void rk4(double ta, double tb, double h, double x, double y, double z, double px, double py,double pz)
{
  std::cout.precision(15);

  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double m1, m2, m3, m4;
  double n1, n2, n3, n4;
  double o1, o2, o3, o4;
  double p1, p2, p3, p4;

  int kk=0;

  const int N = (tb-ta)/h;
  for (int nt = 0; nt < N; ++nt) {
    double t = ta + h*nt;
    
    double k1=0, k2=0, k3=0, k4=0;    
    double l1=0, l2=0, l3=0, l4=0;
    double m1=0, m2=0, m3=0, m4=0;
    double n1=0, n2=0, n3=0, n4=0;
    double o1=0, o2=0, o3=0, o4=0;
    double p1=0, p2=0, p3=0, p4=0;    
    
    k1 = h*f(t, x, y, z, px, py, pz, 0);
    l1 = h*f(t, x, y, z, px, py, pz, 1);
    m1 = h*f(t, x, y, z, px, py, pz, 2);
    n1 = h*f(t, x, y, z, px, py, pz, 3);
    o1 = h*f(t, x, y, z, px, py, pz, 4);
    p1 = h*f(t, x, y, z, px, py, pz, 5);
    
    
    k2= h*f(t + h/2, x+k1/2, y+l1/2, z+m1/2, px+n1/2, py+o1/2, pz+p1/2, 0);
    l2= h*f(t + h/2, x+k1/2, y+l1/2, z+m1/2, px+n1/2, py+o1/2, pz+p1/2, 1);
    m2= h*f(t + h/2, x+k1/2, y+l1/2, z+m1/2, px+n1/2, py+o1/2, pz+p1/2, 2);
    n2= h*f(t + h/2, x+k1/2, y+l1/2, z+m1/2, px+n1/2, py+o1/2, pz+p1/2, 3);
    o2= h*f(t + h/2, x+k1/2, y+l1/2, z+m1/2, px+n1/2, py+o1/2, pz+p1/2, 4);
    p2= h*f(t + h/2, x+k1/2, y+l1/2, z+m1/2, px+n1/2, py+o1/2, pz+p1/2, 3);

    
     
    k3= h*f(t + h/2, x+k2/2, y+l2/2, z+m2/2, px+n2/2, py+o2/2, pz+p2/2, 0);
    l3= h*f(t + h/2, x+k2/2, y+l2/2, z+m2/2, px+n2/2, py+o2/2, pz+p2/2, 1);
    m3= h*f(t + h/2, x+k2/2, y+l2/2, z+m2/2, px+n2/2, py+o2/2, pz+p2/2, 2);
    n3= h*f(t + h/2, x+k2/2, y+l2/2, z+m2/2, px+n2/2, py+o2/2, pz+p2/2, 3);
    o3= h*f(t + h/2, x+k2/2, y+l2/2, z+m2/2, px+n2/2, py+o2/2, pz+p2/2, 4);
    p3= h*f(t + h/2, x+k2/2, y+l2/2, z+m2/2, px+n2/2, py+o2/2, pz+p2/2, 5);

    
    k4= h*f(t + h, x+k3, y+l3, z+m3, px+n3, py+o3, pz+p3, 0);
    l4= h*f(t + h, x+k3, y+l3, z+m3, px+n3, py+o3, pz+p3, 1);
    m4= h*f(t + h, x+k3, y+l3, z+m3, px+n3, py+o3, pz+p3, 2);
    n4= h*f(t + h, x+k3, y+l3, z+m3, px+n3, py+o3, pz+p3, 3);
    o4= h*f(t + h, x+k3, y+l3, z+m3, px+n3, py+o3, pz+p3, 4);
    p4= h*f(t + h, x+k3, y+l3, z+m3, px+n3, py+o3, pz+p3, 5);

    std::cout   << t << "\t\t\t"  << x << "\t" << y <<"\t" << z /* <<"\t\t" << y[3] <<"\t\t" << y[4] <<"\t\t" << y[5] */ << std::endl;

    x = x + (k1 + 2*k2 + 2*k3 + k4)/6.0;
    y = y + (l1 + 2*l2 + 2*l3 + l4)/6.0;
    z = z + (m1 + 2*m2 + 2*m3 + m4)/6.0;
    px = px + (n1 + 2*n2 + 2*n3 + n4)/6.0;
    py = py + (o1 + 2*o2 + 2*o3 + o4)/6.0;
    pz = pz + (p1 + 2*p2 + 2*p3 + p4)/6.0;

    kk=kk+1;
  // std::cout<<kk<<std::endl; 
   //std::cout   << t << "\t\t\t"  << x << "\t" << y <<"\t" << z /* <<"\t\t" << y[3] <<"\t\t" << y[4] <<"\t\t" << y[5] */ << std::endl;
   
  }
}

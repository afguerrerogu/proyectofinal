#include <iostream>
#include <cmath>
#include <vector>

std::vector<double> y ={2.5}; // {x0, y0, z0, px , py, pz} 
std::vector<double> x ={0.0};
std::vector<double> z ={0.0}; 
std::vector<double> px ={0.0}; 
std::vector<double> py ={1.68884}; 
std::vector<double> pz ={0.2}; 

/*std::vector<double> xi ={0.994}; // {x0, y0, u0, vx} 
std::vector<double> yi ={0.0};
std::vector<double> ui ={0.0}; 
std::vector<double> vi ={-2.001585}; 

double const miu= 0.012277471;
double const n=1-miu;*/

const double A= 1.0;
const double C= 1.0;
const double w= 0.25;
const double a= 1.25;
const double b= 1.0;
const double c= 0.75;


double f(double t,double  x,double y,double z,double px, double py,double pz , int id);
//double f(double t, double x, double y, double u, double v, int id);
void rk4(double ta, double tb, double h, std::vector<double> x,std::vector<double> y,std::vector<double> u,std::vector<double> px,std::vector<double> py,std::vector<double> pz);

int main(void){

    rk4(0, 10000 , 0.025, x,y,z,px,py,pz);
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

/*double f(double t, double x, double y, double u, double v, int id){

    if (0==id){
        return u;
    }
    if (1==id){
        return v;
    }
    if (2==id){
        return x+2*v-n*((x+miu)/std::sqrt(pow(pow(x+miu,2)+pow(y,2),3))-miu*((x-n)/std::sqrt(pow(pow(x-n,2)+pow(y,2),3))));
    }
    if (3==id){
        return y-2*u-n*((y)/std::sqrt(pow(pow(x+miu,2)+pow(y,2),3))-miu*((x-n)/std::sqrt(pow(pow(x-n,2)+pow(y,2),3))));
    }
    else
    {
        exit(1);
    }
}*/

void rk4(double ta, double tb, double h, std::vector<double> x,std::vector<double> y,std::vector<double> z,std::vector<double> px ,std::vector<double> py,std::vector<double> pz)
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
    


    k1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt] , 0);
    l1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt] , 1);
    m1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt] , 2);
    n1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt] , 3);
    o1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt] , 4);
    p1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt] , 5);
    
    
    k2= h*f(t + h/2, x[nt]+k1/2, y[nt]+l1/2, z[nt]+m1/2, px[nt]+n1/2, py[nt]+o1/2, pz[nt]+p1/2, 0);
    l2= h*f(t + h/2, x[nt]+k1/2, y[nt]+l1/2, z[nt]+m1/2, px[nt]+n1/2, py[nt]+o1/2, pz[nt]+p1/2, 1);
    m2= h*f(t + h/2, x[nt]+k1/2, y[nt]+l1/2, z[nt]+m1/2, px[nt]+n1/2, py[nt]+o1/2, pz[nt]+p1/2, 2);
    n2= h*f(t + h/2, x[nt]+k1/2, y[nt]+l1/2, z[nt]+m1/2, px[nt]+n1/2, py[nt]+o1/2, pz[nt]+p1/2, 3);
    o2= h*f(t + h/2, x[nt]+k1/2, y[nt]+l1/2, z[nt]+m1/2, px[nt]+n1/2, py[nt]+o1/2, pz[nt]+p1/2, 4);
    p2= h*f(t + h/2, x[nt]+k1/2, y[nt]+l1/2, z[nt]+m1/2, px[nt]+n1/2, py[nt]+o1/2, pz[nt]+p1/2, 5);

    
     
    k3= h*f(t + h/2, x[nt]+k2/2, y[nt]+l2/2, z[nt]+m2/2, px[nt]+n2/2, py[nt]+o2/2, pz[nt]+p2/2, 0);
    l3= h*f(t + h/2, x[nt]+k2/2, y[nt]+l2/2, z[nt]+m2/2, px[nt]+n2/2, py[nt]+o2/2, pz[nt]+p2/2, 1);
    m3= h*f(t + h/2, x[nt]+k2/2, y[nt]+l2/2, z[nt]+m2/2, px[nt]+n2/2, py[nt]+o2/2, pz[nt]+p2/2, 2);
    n3= h*f(t + h/2, x[nt]+k2/2, y[nt]+l2/2, z[nt]+m2/2, px[nt]+n2/2, py[nt]+o2/2, pz[nt]+p2/2, 3);
    o3= h*f(t + h/2, x[nt]+k2/2, y[nt]+l2/2, z[nt]+m2/2, px[nt]+n2/2, py[nt]+o2/2, pz[nt]+p2/2, 4);
    p3= h*f(t + h/2, x[nt]+k2/2, y[nt]+l2/2, z[nt]+m2/2, px[nt]+n2/2, py[nt]+o2/2, pz[nt]+p2/2, 5);
    
    k4= h*f(t + h, x[nt]+k3, y[nt]+l3, z[nt]+m3, px[nt]+n3, py[nt]+o3, pz[nt]+p3, 0);
    l4= h*f(t + h, x[nt]+k3, y[nt]+l3, z[nt]+m3, px[nt]+n3, py[nt]+o3, pz[nt]+p3, 1);
    m4= h*f(t + h, x[nt]+k3, y[nt]+l3, z[nt]+m3, px[nt]+n3, py[nt]+o3, pz[nt]+p3, 2);
    n4= h*f(t + h, x[nt]+k3, y[nt]+l3, z[nt]+m3, px[nt]+n3, py[nt]+o3, pz[nt]+p3, 3);
    o4= h*f(t + h, x[nt]+k3, y[nt]+l3, z[nt]+m3, px[nt]+n3, py[nt]+o3, pz[nt]+p3, 4);
    p4= h*f(t + h, x[nt]+k3, y[nt]+l3, z[nt]+m3, px[nt]+n3, py[nt]+o3, pz[nt]+p3, 5);

 
    x[nt+1] = x[nt] + (k1 + 2*k2 + 2*k3 + k4)/6.0;
    y[nt+1] = y[nt] + (l1 + 2*l2 + 2*l3 + l4)/6.0;
    z[nt+1] = z[nt] + (m1 + 2*m2 + 2*m3 + m4)/6.0;
    px[nt+1] = px[nt] + (n1 + 2*n2 + 2*n3 + n4)/6.0;
    py[nt+1] = py[nt] + (o1 + 2*o2 + 2*o3 + o4)/6.0;
    pz[nt+1] = pz[nt] + (p1 + 2*p2 + 2*p3 + p4)/6.0;

    if (f(nt, x[nt+1], y[nt+1],z[nt+1], px[nt+1],py[nt+1],pz[nt+1],1)>0 && y[nt+1]>0)
    {
     std::cout  /* << t << "\t\t\t" */ << x[nt+1] << "\t" << y[nt+1] <<"\t" << z[nt+1] /* <<"\t\t" << y[3] <<"\t\t" << y[4] <<"\t\t" << y[5] */ << std::endl;
    }


    //kk=kk+1;
    //std::cout<<kk<<std::endl;
  
   
  }
}
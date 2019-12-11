#include <iostream>
#include <cmath>
#include <vector>

std::vector<double> y = {2.5, 0 , 0 , 0 , 1.68884, 0.2 }; // {x0, y0, z0, px , py, pz} 
//std::vector<double> y = {0.994, 0 , 0 ,-2.001585 }; // {x0, y0, u0, vx }

const double A= 1.0;
const double C= 1.0;
const double w= 0.25;
const double a= 1.25;
const double b= 1.0;
const double c= 0.75;

/*double const miu= 0.012277471;
double const n=1-miu;*/

double f(double t, const std::vector<double>  y, int id);
void rk4(double ta, double tb, double h, std::vector<double> & y);

int main(void){

    rk4(0, 10000 , 0.025, y);
    return 0;
}

double f(double t, const std::vector<double>  y, int id){

    if (0==id){
        return y[3]+w*y[1];
    }
    if (1==id){
        return y[4]-w*y[0];
    }
    if (2==id){
        return y[5]; 
    }
    if (3==id){
        return w*y[4]-(2*A*y[0])/(pow(a,2)*(C+pow(y[0]/a ,2)+pow(y[1]/b ,2)+pow(y[2]/c ,2)));
    }
    if (4==id){
        return -w*y[3]-(2*A*y[1])/(pow(b,2)*(C+pow(y[0]/a ,2)+pow(y[1]/b ,2)+pow(y[2]/c ,2)));
    }
    if (5==id){
        return -(2*A*y[2])/(pow(c,2)*(C+pow(y[0]/a ,2)+pow(y[1]/b ,2)+pow(y[2]/c ,2)));
    }
    else
    {
        exit(1);
    }
}

/*double f(double t, const std::vector<double>  y, int id){

    if (0==id){
        return y[2];
    }
    if (1==id){
        return y[3];
    }
    if (2==id){
        return y[0]+2*y[3]-n*((y[0]+miu)/std::sqrt(pow(pow(y[0]+miu,2)+pow(y[1],2),3))-miu*((y[0]-n)/std::sqrt(pow(pow(y[0]-n,2)+pow(y[1],2),3))));
    }
    if (3==id){
        return y[1]-2*y[2]-n*((y[1])/std::sqrt(pow(pow(y[0]+miu,2)+pow(y[1],2),3))-miu*((y[0]-n)/std::sqrt(pow(pow(y[0]-n,2)+pow(y[1],2),3))));
    }
    else
    {
        exit(1);
    }
}*/

void rk4(double ta, double tb, double h, std::vector<double> & y)
{
  std::cout.precision(15);

  std::vector<double> k1, k2, k3, k4, aux;
  k1.resize(y.size());
  k2.resize(y.size());
  k3.resize(y.size());
  k4.resize(y.size());
  aux.resize(y.size());

  const int N = (tb-ta)/h;
  for (int nt = 0; nt < N; ++nt) {
    double t = ta + h*nt;
    
    // k1
    for(int ii = 0; ii < y.size(); ++ii) {
      k1[ii] = h*f(t, y, ii);
    }

    
    // k2 aux
    for(int ii = 0; ii < y.size(); ++ii) {
      aux[ii] = y[ii] + k1[ii]/2;
    }
    //k2
    for(int ii = 0; ii < y.size(); ++ii) {
      k2[ii] = h*f(t + h/2, aux, ii);
    }
    

    // k3 aux
    for(int ii = 0; ii < y.size(); ++ii) {
      aux[ii] = y[ii] + k2[ii]/2;
    }
    //k3
    for(int ii = 0; ii < y.size(); ++ii) {
      k3[ii] = h*f(t + h/2, aux, ii);
    }
    

    // k4 aux
    for(int ii = 0; ii < y.size(); ++ii) {
      aux[ii] = y[ii] + k3[ii];
    }
    //k4
    for(int ii = 0; ii < y.size(); ++ii) {
      k4[ii] = h*f(t + h, aux, ii);
    } 

    // write new y
    for(int ii = 0; ii < y.size(); ++ii) {
      y[ii] = y[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
    }

    if (f(nt, y, 1)>0 && y[1]> 0)
    {
      std::cout  /* << t << "\t\t\t" */ << y[0] << "\t" << y[1] <<"\t" << y[2] /* <<"\t\t" << y[3] <<"\t\t" << y[4] <<"\t\t" << y[5] */ << std::endl;
    }
    
   
   
  }
}

import numpy as np
import matplotlib.pyplot as plt

ti=0;
tf=10000;
h=0.1;
N=((tf-ti)/h);

#CONDICIONES INICIALES
x0 = 0;
y0 = 2.5;
z0 = 0;
px0 = 0;
py0 = 1.68884;
pz0 = 0.2;

A= 1.0;
C= 1.0;
w= 0.25;
a= 1.25;
b= 1.0;
c= 0.75;

def f(t, x, y, z, px, py, pz , id):
    if (0==id):
        return px+w*y
    if (1==id):
        return py-w*x
    if (2==id):
        return pz 
    if (3==id):
        return w*py-(2*A*x)/(pow(a,2)*(C+pow(x/a ,2)+pow(y/b ,2)+pow(z/c ,2)))
    if (4==id):
        return -w*px-(2*A*y)/(pow(b,2)*(C+pow(x/a ,2)+pow(y/b ,2)+pow(z/c ,2)))
    if (5==id):
        return -(2*A*z)/(pow(c,2)*(C+pow(x/a ,2)+pow(y/b ,2)+pow(z/c ,2)))
    else:
        exit(1)




def runge_kutta_sis(f,x,y,z, px, py, pz, ti,tf,h):
    t = np.arange(ti,tf+h,h)
    n = len(t)
    x = np.zeros(10000000);
    y = np.zeros(10000000);
    z = np.zeros(10000000);
    px =np.zeros(10000000);
    py =np.zeros(10000000);
    pz =np.zeros(10000000);
    
    x[0] = x0; y[0] = y0; z[0] = z0; px[0] = px0; py[0] = py0; pz[0] = pz0;
    
    for nt in range(n-1):

        k1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt], 0);
        l1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt], 1);
        m1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt], 2);
        n1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt], 3);
        o1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt], 4);
        p1 = h*f(t, x[nt], y[nt], z[nt], px[nt], py[nt], pz[nt], 5);
    
    
        k2= h*f(t + h/2, x[nt]+k1/2, y[nt]+l1/2, z[nt]+m1/2, px[nt]+n1/2, py[nt]+o1/2, pz[nt]+p1/2, 0);
        l2=h*f(t + h/2, x[nt]+k1/2, y[nt]+l1/2, z[nt]+m1/2, px[nt]+n1/2, py[nt]+o1/2, pz[nt]+p1/2, 1);
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

    
        print (x[nt+1],y[nt+1],z[nt+1])

runge_kutta_sis(f, x0,y0,z0,px0,py0,pz0, ti,tf,h)

#fixed end beams

import numpy, matplotlib.pyplot

kb=0.5 #usualy the bending stiffness of PP (plastic) ranges from 1 to 10 GPa. 1 GPa=1000N/mm^2=0.1N/cm^2
ks=1.4
steps=100
h=1./steps
a0=6*kb+2*ks*h**2
a1=-4*kb-ks*h**2
a2=kb
length=30. #unit:cm


def Fixed_End_Beam(delta):
    A=numpy.zeros((steps+1,steps+1))
    b=numpy.zeros(steps+1)

    A[0][0]=1.
    b[0]=0.

    A[1][0]=a1
    A[1][1]=a0+a2
    A[1][2]=a1
    A[1][3]=a2

    N = int(steps / 2)
    for i in range(2,N):
        A[i][i-2]=a2
        A[i][i-1]=a1
        A[i][i]=a0
        A[i][i+1]=a1
        A[i][i+2]=a2
        b[i]=0

    A[N][N]=1.
    b[N]=delta

    for i in range(N+1,steps-1):
        A[i][i-2]=a2
        A[i][i-1]=a1
        A[i][i]=a0
        A[i][i+1]=a1
        A[i][i+2]=a2
        b[i]=0

    A[steps - 1][steps - 3] = a2
    A[steps - 1][steps - 2] = a1
    A[steps - 1][steps - 1] = a0 + a2
    A[steps-1][steps]=a1
    b[steps-1]=0

    A[steps][steps]=1.
    b[steps]=0
    x=numpy.linalg.solve(A,b)
    return x


def Cantilever_Beam(delta,f):
    A=numpy.zeros((steps+1,steps+1))
    b=numpy.zeros(steps+1)

    A[0][0]=1.
    b[0]=0.

    A[1][0]=a1
    A[1][1]=a0+a2
    A[1][2]=a1
    A[1][3]=a2

    for i in range(2,steps-2):
        A[i][i-2]=a2
        A[i][i-1]=a1
        A[i][i]=a0
        A[i][i+1]=a1
        A[i][i+2]=a2
        b[i]=0

    A[steps - 2][steps - 4] = a2
    A[steps - 2][steps - 3] = a1
    A[steps - 2][steps - 2] = a0 - a2
    A[steps-2][steps-1]=a1+2*a2
    b[steps-2]=0

    A[steps-1][steps-3]=-1.
    A[steps-1][steps-2]=3.
    A[steps-1][steps-1]=-3.
    A[steps-1][steps]=1
    b[steps-1]=h**3*f/kb

    A[steps][steps]=1.
    b[steps]=delta
    x=numpy.linalg.solve(A,b)
    return x


def Simply_Supported_Beam(delta):
    A=numpy.zeros((steps+1,steps+1))
    b=numpy.zeros(steps+1)

    A[0][0]=1.
    b[0]=0.

    A[1][1]=a0-a1
    A[1][2]=a2
    A[1][3]=a1
    b[1]=0.

    N=int(steps / 2)
    for i in range(2,N-2):
        A[i][i-2]=a2
        A[i][i-1]=a1
        A[i][i]=a0
        A[i][i+1]=a1
        A[i][i+2]=a2
        b[i]=0

    A[N- 2][N - 2] = a2
    A[N - 2][N - 1]=a1
    A[N - 2][N]=a0
    A[N - 2][N + 1] = a1
    b[N-2]=-a2*delta

    A[N-1][N-3]=a2
    A[N-1][N-2]=a1
    A[N-1][N-1]=a0
    A[N-1][N+1]=a2
    b[N-1]=-a1*delta

    A[N][N]=1.
    b[N]=delta

    A[N+1][N-1]=a2
    A[N+1][N+1]=a0
    A[N+1][N+2]=a1
    A[N + 1][N + 3] = a2
    b[N+1]=-a1*delta

    A[N+2][N+1]=a1
    A[N+2][N+2]=a0
    A[N+2][N+3]=a1
    A[N+2][N+4]=a2
    b[N+2]=-a2*delta

    for i in range(N+3,steps-1):
        A[i][i-2]=a2
        A[i][i-1]=a1
        A[i][i]=a0
        A[i][i+1]=a1
        A[i][i+2]=a2
        b[i]=0

    A[steps - 1][steps - 3] = a2
    A[steps - 1][steps - 2] = a1
    A[steps-1][steps-1]=a0-a2
    b[steps-1]=0

    A[steps][steps]=1.
    b[steps]=0.

    x=numpy.linalg.solve(A,b)
    return x


def plot_simply_supported_beam(x,y):
    matplotlib.pyplot.axis([0, length, -20, 5])
    matplotlib.pyplot.plot(x,y)
    matplotlib.pyplot.xlabel('Length in cm')
    matplotlib.pyplot.ylabel('Height in cm')
    matplotlib.pyplot.title('Simply Supported Beam')
    matplotlib.pyplot.axis('equal')
    matplotlib.pyplot.text(x=0,y=2,s='delta='+str(delta))
    matplotlib.pyplot.show()


def plot_fixed_end_beam(x,y):
    matplotlib.pyplot.axis([0, length, -20, 5])
    matplotlib.pyplot.plot(x,y)
    matplotlib.pyplot.xlabel('Length in cm')
    matplotlib.pyplot.ylabel('Height in cm')
    matplotlib.pyplot.title('Fixed End Beam')
    matplotlib.pyplot.axis('equal')
    matplotlib.pyplot.text(x=0,y=2,s='delta='+str(delta))
    matplotlib.pyplot.show()


def plot_cantilever_beam(x, y):
    matplotlib.pyplot.axis([0, length, -20, 5])
    matplotlib.pyplot.plot(x, y)
    matplotlib.pyplot.xlabel('Length in cm')
    matplotlib.pyplot.ylabel('Height in cm')
    matplotlib.pyplot.title('Cantilever Beam')
    matplotlib.pyplot.axis('equal')
    matplotlib.pyplot.text(x=0, y=2, s='delta=' + str(delta))
    matplotlib.pyplot.show()


if __name__ == "__main__":
    x=numpy.zeros(steps+1)

    delta = -1.06
    for i in range(steps+1):
        x[i]=length*i/steps
    y=Simply_Supported_Beam(delta)
    plot_simply_supported_beam(x,y)

    delta = -1.5
    length=40.8
    for i in range(steps+1):
        x[i]=length*i/steps
    y = Fixed_End_Beam(delta)
    plot_fixed_end_beam(x,y)

    delta = -7.2
    length=39
    f=40./1000.*9.8
    for i in range(steps+1):
        x[i]=length*i/steps
    y = Cantilever_Beam(delta,f)
    plot_cantilever_beam(x,y)
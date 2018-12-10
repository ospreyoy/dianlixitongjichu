#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
using namespace std;
int node_num;

double qian_P(double P2out,double Q2out,double U2,double U1,double r0,double x0,double b0,double l)
{
    double R;
    double X;
    double B;
    if(l==0)
    {
        R=r0;
        X=x0;
        B=b0;
    }
    else
    {
        R=r0*l*(1.0-((x0*b0*l*l)/3.0));
        X=x0*l*(1-((l*l*((x0*b0)-((r0*r0*b0)/x0)))/6));
        B=b0*l*(1+((x0*b0*l*l)/12));
    }
    double P2;
    double Q2;
    P2=P2out;
    Q2=Q2out-((U2*U2*B)/4);

    double dP2;
    double dQ2;
    dP2=R*((P2*P2+Q2*Q2)/(U2*U2));
    dQ2=X*((P2*P2+Q2*Q2)/(U2*U2));

    double P1;
    double Q1;
    P1=P2+dP2;
    Q1=Q2+dQ2;

    double P1in;
    double Q1in;
    P1in=P1;
    Q1in=Q1-((U1*U1*B)/4);
    return P1in;
}
double qian_Q(double P2out,double Q2out,double U2,double U1,double r0,double x0,double b0,double l)
{
    double R;
    double X;
    double B;

    if(l==0)
    {
        R=r0;
        X=x0;
        B=b0;
    }
    else
    {
        R=r0*l*(1.0-((x0*b0*l*l)/3.0));
        X=x0*l*(1-((l*l*((x0*b0)-((r0*r0*b0)/x0)))/6));
        B=b0*l*(1+((x0*b0*l*l)/12));
    }
    double P2;
    double Q2;
    P2=P2out;
    Q2=Q2out-((U2*U2*B)/4);

    double dP2;
    double dQ2;
    dP2=R*((P2*P2+Q2*Q2)/(U2*U2));
    dQ2=X*((P2*P2+Q2*Q2)/(U2*U2));

    double P1;
    double Q1;
    P1=P2+dP2;
    Q1=Q2+dQ2;

    double P1in;
    double Q1in;
    P1in=P1;
    Q1in=Q1-((U1*U1*B)/4);
    return Q1in;
}
double hou_U(double P1in,double Q1in,double U1,double r0,double x0,double b0,double l)
{
    double R;
    double X;
    double B;
    if(l==0)
    {
        R=r0;
        X=x0;
        B=b0;
    }
    else
    {
        R=r0*l*(1.0-((x0*b0*l*l)/3.0));
        X=x0*l*(1-((l*l*((x0*b0)-((r0*r0*b0)/x0)))/6));
        B=b0*l*(1+((x0*b0*l*l)/12));
    }
    double P1;
    double Q1;
    P1=P1in;
    Q1=Q1in+((U1*U1*B)/4);

    double DU1;
    double dU1;
    DU1=(P1*R+Q1*X)/U1;
    dU1=(P1*X-Q1*R)/U1;

    double U2;
    U2=sqrt((U1-DU1)*(U1-DU1)+dU1*dU1);
    return U2;
}
int main()
{
    ifstream in("data.txt");
    double U1;//首端真实电压
    double UN;//线路额定电压

    in>>node_num;
    in>>U1;
    in>>UN;

    double *r0=new double[node_num];//定义各段参数 电阻
    double *x0=new double[node_num];//电抗
    double *b0=new double[node_num];//电纳
    double *l=new double[node_num]; //定义各个节点之间长度 如为集总参数则置0

    for(int xunhuan0=0;xunhuan0<node_num;xunhuan0++) //参数读入赋值
    {
        in>>r0[xunhuan0];
        in>>x0[xunhuan0];
        in>>b0[xunhuan0];
        in>>l[xunhuan0];
    }

    double *Pout=new double[node_num];//定义各节点输出功率
    double *Qout=new double[node_num];
    for(int xunhuan2=0;xunhuan2<node_num;xunhuan2++)
    {
        in>>Pout[xunhuan2];
        in>>Qout[xunhuan2];
    }
    int diedai_num;
    in>>diedai_num;

    double *Pchuan=new double[node_num+1];
    double *Qchuan=new double[node_num+1];
    Pchuan[node_num]=0;
    Qchuan[node_num]=0;
    double *U=new double[node_num+1];
    U[0]=U1;
    for(int i=0;i<node_num;i++)
    {
        U[i+1]=UN;
    }


    ofstream out("out.txt");
    for(int k=0;k<diedai_num;k++)   //迭代循环
    {
               //一次迭代
    for(int i=0;i<node_num;i++)
    {
        Pchuan[node_num-i-1]=qian_P(Pout[node_num-i-1]+Pchuan[node_num-i],Qout[node_num-i-1]+Qchuan[node_num-i],U[node_num-i],U[node_num-i-1],r0[node_num-i-1],x0[node_num-i-1],b0[node_num-i-1],l[node_num-i-1]);
        Qchuan[node_num-i-1]=qian_Q(Pout[node_num-i-1]+Pchuan[node_num-i],Qout[node_num-i-1]+Qchuan[node_num-i],U[node_num-i],U[node_num-i-1],r0[node_num-i-1],x0[node_num-i-1],b0[node_num-i-1],l[node_num-i-1]);
    }
    for(int i=0;i<node_num;i++)
    {
        U[i+1]=hou_U(Pchuan[i],Qchuan[i],U[i],r0[i],x0[i],b0[i],l[i]);
    }
    out<<"第"<<k+1<<"次迭代各节点电压为：";
    for(int i=0;i<node_num;i++)
    {
        out<<"U("<<i+1<<")="<<U[i+1]<<"    ";
    }
    out<<endl;
    out<<"各节点向下传送功率为："<<"  ";
    for(int i=0;i<node_num;i++)
    {
        out<<"S("<<i<<"→"<<i+1<<")="<<Pchuan[i]<<"+j"<<Qchuan[i]<<"    ";
    }
    out<<endl<<endl;
    }
    delete []r0;
    delete []x0;
    delete []b0;
    delete []l;
    delete []Pout;
    delete []Qout;
    delete []Pchuan;
    delete []Qchuan;
    delete []U;
    return 0;
}

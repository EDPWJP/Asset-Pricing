#include<iostream>
#include<cmath>
#include<vector>
#include<cnormdist>
using namespace std;

double max(const double &a, const double &b)
{
	if(x>y) return x;
	else return y;
}

double bond_price_discrete(const vector<double> &times, const vector<double> &cashflows, const double &r)//Bond Price
{
	double price=0.0;
	for(int i=0;i<times.size();i++)
		{
			price+=cashflows[i]/pow(1.0+r,times[i]);
		}
	return price;
}

double bonds_yield_to_maturity_discrete(const vector<double> &times, const vector<double> &cashflows, const double &bondprice)//Yield to Maturity
{
	double bot=0.0, top=1.0,y=0.5;
	double ACCURACY=1e-5;//least tolerated diff
	while(bond_price_discrete(times,cashflows,top)>bondprice)
		{top=top*2;}
	while(fabs(bond_price_discrete(times,cashflows,y)-bondprice)>ACCURACY)
	{
		if (bond_price_discrete(times,cashflows,y)>bondprice)
			{
				bot=y;//update bot to y if y is too small
				y=(top+y)/2;//update y to narrow down
			}
		else 
			{
				top=y;//update top to y if y is too large
				y=(bot+y)/2;//update y to narrow down
			}
	}
	return y;
}

double average(const vector<double> &times, const vector<double> &cashflows)//Historical Average
{
	double sum=0.0; 
	vector<double> yield(times.size());
	for(int i=1;i<times.size();i++)
		{
			yield[i]=(cashflows[i]-cashflows[i-1])/cashflows[i-1];
			sum+=yield[i];
		}
	double average=sum/(times.size()-1);
	return average;
}

double standard_deviation(const vector<double> &times, const vector<double> &cashflows)//Historical std
{
	double sum=0.0;
	vector<double> yield(times.size());
	double aver=average(times,cashflows);
	for(int i=1;i<times.size();i++)
		{
			yield[i]=(cashflows[i]-cashflows[i-1])/cashflows[i-1];
			sum+=pow(yield[i]-aver,2);
		}
	double sd=pow(sum/(times.size()-2),0.5);//sample var use n-1
	return sd;
}

double N(const double &x)//N(x) Normal CDF
{
	if(x>6.0){return 1.0;}
	if(x<-6.0){return 0.0;}
	double b1=0.31938153;
	double b2=-0.356563782;
	double b3=1.781477937;
	double b4=-1.821255978;
	double b5=1.330274429;
	double p=0.2316419;
	double c2=0.3989423;
	double a=fabs(x);
	double t=1.0/(1.0+a*p);
	double n=1.0-c2*exp(-x*x/2)*t*((((b5*t+b4)*t+b3)*t+b2)*t+b1);//use approximation of N(x)
	if(x<0.0){n=1.0-n;}
	return n;
}

double option_price_call_black_scholes(const double &S, const double &X, const double &r, const double &sigma, const double &time)//BSM call
{
	double d1=(log(S/X)+r*time)/sigma*sqrt(time)+0.5*sigma*sqrt(time);
	double d2=d1-sigma*sqrt(time);
	double c=S*N(d1)-X*exp(-r*time)*N(d2);//Black-Scholes Model
	return c;
}

double option_price_put_black_scholes(const double &S, const double &X, const double &r, const double &sigma, const double &time)//BSM put
{
	double d1=(log(S/X)+r*time)/sigma*sqrt(time)+0.5*sigma*sqrt(time);
	double d2=d1-sigma*sqrt(time);
	double p=X*exp(-r*time)*N(-d2)-S*N(-d1);//Black-Scholes Model
	return p;
}

void option_greek_letters(const double &S, const double &X, const double &r, const double &sigma, const double &time)//calculate greek letters of BSM call option
{
	double d1=(log(S/X)+r*time)/sigma*sqrt(time)+0.5*sigma*sqrt(time);
	double d2=d1-sigma*sqrt(time);
	Delta=N(d1);
	Gamma=n(d1)/(S*sigma*sqrt(time));
	Theta=-(S*sigma*n(d1))/(2*sqrt(time))-r*X*exp(-r*time)*N(d2);
	Vega=S*sqrt(time)*n(d1);
	Rho=X*time*exp(-r*time)*N(d2);
}

double option_price_call_monte_carlo(const double &S, const double &X, const double &r, const double &sigma, const double &time, const double &no_sims)//may use control variate or antithesis to optimize
{
	double R=(r-0.5*pow(sigma,2))*time; 
	double SD=sigma*sqrt(time);
	double sum_payoff=0.0;
	for(int i=1;i<=no_sims;i++)
	{
		double S_T=S*exp(R+SD*random_normal());// explicit solution using simulation path
		sum_payoff+=max(0.0,S_T-X);
	}
	call_price=exp(-r*time)*(sum_payoff/no_sims);	
}

double option_price_put_monte_carlo(const double &S, const double &X, const double &r, const double &sigma, const double &time, const double &no_sims)
{
	double R=(r-0.5*pow(sigma,2))*time;
	double SD=sigma*sqrt(time);
	double sum_payoff=0.0;
	for(int i=1;i<=no_sims;i++)
	{
		double S_T=S*exp(R+SD*random_normal());
		sum_payoff+=max(0.0,X-S_T);
	}
	put_price=exp(-r*time)*(sum_payoff/no_sims);	
}

double option_price_call_european_binomial(const double &S, const double &X, const double &r, const double &sigma, const double &time, const double &steps)
{
	double R=exp(r*time/steps);
 	double Rinv=1/R;
	double u=exp(sigma*sqrt(time/steps));
 	double d=1/u;
 	double p_up=(R-d)/(u-d);
 	double p_down=1-p_up;
 	vector<double> stock_price(steps+1);
 	vector<double> call_value(steps+1);
 	stock_price[0]=S*pow(d,steps);
 	for (int i=1;i<=steps;i++) 
	 {
	 	stock_price[i]=u*u*stock_price[i-1];
	 }
	 for(int i =0,i<=steps,i++)
	 {
	 	call_value[i]=max(0.0,stock_price[i]-X);
	 }
	 //initialize the last column
	 for(int step=steps-1,step>=0,step--)//loop through each step from the back
	 {
	 	for(int i=0,i<=step,i++)//loop through every stock price possibility in that step
		 {
		 	stock_price[i]=d*stock_price[i+1];
		 	call_value[i]=Rinv*(p_up*call_value[i+1]+p_down*call_value[i]);
		 }
	 }
	 return call_value[0];
}

double option_price_put_european_binomial(const double &S, const double &X, const double &r, const double &sigma, const double &time, const double &steps)
{
	double R=exp(r*time/steps);
 	double Rinv=1/R;
	double u=exp(sigma*sqrt(time/steps));
 	double d=1/u;
 	double p_up=(R-d)/(u-d);
 	double p_down=1-p_up;
 	vector<double> stock_price(steps+1);
 	vector<double> put_value(steps+1);
 	stock_price[0]=S*pow(d,steps);
 	for (int i=1;i<=steps;i++) 
	 {
	 	stock_price[i]=u*u*stock_price[i-1];
	 }
	 for(int i =0;i<=steps;i++)
	 {
	 	put_value[i]=max(0.0,X-stock_price[i]);
	 }
	 for(int step=steps-1;step>=0;step--)//loop through each step
	 {
	 	for(int i=0;i<=step;i++)//loop through every stock price possibility in that step
		 {
		 	stock_price[i]=d*stock_price[i+1];
		 	put_value[i]=Rinv*(p_up*put_value[i+1]+p_down*put_value[i]);
		 }
	 }
	 return put_value[0];
}

double option_price_call_american_binomial(const double &S, const double &X, const double &r, const double &sigma, const double &time, const double &steps)
{
	double R=exp(r*time/steps);
 	double Rinv=1/R;
	double u=exp(sigma*sqrt(time/steps));
 	double d=1/u;
 	double p_up=(R-d)/(u-d);
 	double p_down=1-p_up;
 	vector<double> stock_price(steps+1);
 	vector<double> call_value(steps+1);
 	stock_price[0]=S*pow(d,steps);
 	for (int i=1;i<=steps;i++) 
	 {
	 	stock_price[i]=u*u*stock_price[i-1];
	 }
	 for(int i =0,i<=steps,i++)
	 {
	 	call_value[i]=max(0.0,stock_price[i]-X);
	 }
	 for(int step=steps-1,step>=0,step--)//loop through each step
	 {
	 	for(int i=0,i<=step,i++)//loop through every stock price possibility in that step
		 {
		 	stock_price[i]=d*stock_price[i+1];
		 	call_value[i]=max(stock_price[i]-X,Rinv*(p_up*call_value[i+1]+p_down*call_value[i]));//check early exercise
		 }
	 }
	 return call_value[0];
}

double option_price_put_american_binomial(const double &S, const double &X, const double &r, const double &sigma, const double &time, const double &steps)
{
	double R=exp(r*time/steps);
 	double Rinv=1/R;
	double u=exp(sigma*sqrt(time/steps));
 	double d=1/u;
 	double p_up=(R-d)/(u-d);
 	double p_down=1-p_up;
 	vector<double> stock_price(steps+1);
 	vector<double> put_value(steps+1);
 	stock_price[0]=S*pow(d,steps);
 	for (int i=1;i<=steps;i++) 
	 {
	 	stock_price[i]=u*u*stock_price[i-1];
	 }
	 for(int i =0,i<=steps,i++)
	 {
	 	put_value[i]=max(0.0,X-stock_price[i]);
	 }
	 for(int step=steps-1,step>=0,step--)//loop through each step
	 {
	 	for(int i=0,i<=step,i++)//loop through every stock price possibility in that step
		 {
		 	stock_price[i]=d*stock_price[i+1];
		 	put_value[i]=max(X-stock_price[i],Rinv*(p_up*put_value[i+1]+p_down*put_value[i]));//check early exercise
		 }
	 }
	 return put_value[0];
}

double option_price_call_european_finite_diff_explicit(const double &S, const double &X, const double &r, const double &sigma, const double &time, const double &no_S_steps, const double &no_t_steps)
{
	int N=no_t_steps;
	double delta_t=time/N;
	int M;
	if(no_S_steps%2==0){M=no_S_steps;}
	else{M=no_S_steps+1;}
	double delta_S=2.0*S/M;
	vector<double> S_values(M+1);
	for(int i=0;i<=M;i++)
	{
		S_values[i]=i*delta_S;
	}
	//set steps delta and values for t and S
	vector <double> a(M);
	vector <double> b(M);
	vector <double> c(M);
	for(int i=0;i<M;i++)
	{
		a[i]=0.5*i*(-r+sigma*sigma*i)*delta_t/(1.0+r*delta_t);
		b[i]=(1.0-sigma*sigma*i*i*delta_t)/(1.0+r*delta_t);
		c[i]=0.5*i*(r+sigma*sigma*i)*delta_t/(1.0+r*delta_t);
	}
	//set parameters a b c for backward difference equation
	vector <double> f_next(M+1);
	for(int i=0;i<=M;i++)
	{
		f_next[i]=max(0.0,S_values[i]-X);
	}
	//initialize last row
	vector <double> f(M)
	for(int t=N-1;t>=0,t--)
	{
		f[0]=0; //assign 0 to edge points to maintain constant size of column
		for(int m=1,m<M,m++)
		{
			f[m]=a[m]*f_next[m-1]+b[m]*f_next[m]+c[m]*f_next[m+1];//use f_next to obtain f
		}
		f[M]=0; //assign 0 to edge points
		for(int i=0;i<=M;i++)
		{
			f_next[i]=f[i];//update next row
		}
	}
return f[M/2];
}

double option_price_call_american_finite_diff_explicit(const double &S, const double &X, const double &r, const double &sigma, const double &time, const double &no_S_steps, const double &no_t_steps)
{
	int N=no_t_steps;
	double delta_t=time/N;
	int M;
	if(no_S_steps%2==0){M=no_S_steps;}
	else{M=no_S_steps+1;}
	double delta_S=2.0*S/M;
	vector<double> S_values(M+1);
	for(int i=0;i<=M;i++)
	{
		S_values[i]=i*delta_S;
	}
	//set steps delta and values for t and S
	vector <double> a(M);
	vector <double> b(M);
	vector <double> c(M);
	for(int i=0;i<M;i++)
	{
		a[i]=0.5*i*(-r+sigma*sigma*i)*delta_t/(1.0+r*delta_t);
		b[i]=(1.0-sigma*sigma*i*i*delta_t)/(1.0+r*delta_t);
		c[i]=0.5*i*(r+sigma*sigma*i)*delta_t/(1.0+r*delta_t);
	}
	//set parameters a b c for backward difference equation
	vector <double> f_next(M+1);
	for(int i=0;i<=M;i++)
	{
		f_next[i]=max(0.0,S_values[i]-X);
	}
	vector <double> f(M)
	for(int t=N-1;t>=0,t--)
	{
		f[0]=0;
		for(int m=1,m<M,m++)
		{
			f[m]=max(a[m]*f_next[m-1]+b[m]*f_next[m]+c[m]*f_next[m+1],S_values[m]-X);//check early exercise
		}
		f[M]=0;
		for(int i=0;i<=M;i++)
		{
			f_next[i]=f[i];//update next row
		}
	}
return f[M/2];
}

main()
{	double bondprice=256.0, r=0.1;
	vector<double>times;
	times.push_back(1);
	times.push_back(2);
	times.push_back(3);
	times.push_back(4);
	vector<double>cashflows;
	cashflows.push_back(100.0);
	cashflows.push_back(200.0);
	cashflows.push_back(300.0);
	cashflows.push_back(400.0);
	cout<<bond_price_discrete(times,cashflows,r)<<endl;
	cout<<bonds_yield_to_maturity_discrete(times,cashflows,bondprice)<<endl;
	cout<<average(times,cashflows)<<endl;
	cout<<standard_deviation(times,cashflows)<<endl;	
}

/*
Implementation of the class WaveEq1D for solving the 1D time-dependent scalar wave equation

Nguyen Trung Thanh, UNCC 2014
*/

#include "include/WaveEq1D.h"

using namespace std;

//implementation of member functions: 

WaveEq1D::WaveEq1D() // default constructor
{ 
	m_grid = Grid1D();
	m_MaxTime = 0.0;
	m_TimeIncWaveExcited = 0.0;
	m_dt = 0.0;
	m_NrTimeSteps = 0;
	m_Nnodes = 0;
	Vec_real x(1); x = 0.0;  
	m_h = x; 
	m_Coef = x; 
	m_IncWaveFunc = x; 
	m_MassMatrix_l = x;
	m_MassMatrix_m = x;
	m_MassMatrix_u = x; 
	m_StiffMatrix_l = x;
	m_StiffMatrix_m = x;
	m_StiffMatrix_u = x;
	m_Solution_previous = x;
	m_Solution_current = x;
	m_Solution_next = x;
	m_Solution.newsize(1,1); m_Solution = 0.0; 


}
WaveEq1D::WaveEq1D(Grid1D grid, double MaxTime, double TimeIncWaveExcited, double dt)
{
	m_grid = grid;
	m_MaxTime = MaxTime;
	m_TimeIncWaveExcited = TimeIncWaveExcited;
	m_dt = dt;
	m_NrTimeSteps = round(MaxTime/dt) + 1;
	m_Nnodes = grid.get_nnodes();
	Vec_real x(1); x = 0.0;  
	Vec_real h(m_Nnodes-1); 
	for (int i = 0; i < m_Nnodes-1; i++)
		h(i) = m_grid.get_meshsize(i);
	m_h = h; 
	m_Coef = x; 
	m_IncWaveFunc = x; 
	m_MassMatrix_l = x;
	m_MassMatrix_m = x;
	m_MassMatrix_u = x; 
	m_StiffMatrix_l = x;
	m_StiffMatrix_m = x;
	m_StiffMatrix_u = x;
	m_Solution_previous = x;
	m_Solution_current = x;
	m_Solution_next = x;
	m_Solution.newsize(1,1); m_Solution = 0.0; 
	
}

WaveEq1D::WaveEq1D(Grid1D grid, double MaxTime, double TimeIncWaveExcited, double dt, Vec_real Coefficient, Vec_real IncWaveFunc)
{
	m_grid = grid;
	m_MaxTime = MaxTime;
	m_TimeIncWaveExcited = TimeIncWaveExcited;
	m_dt = dt;
	m_NrTimeSteps = round(MaxTime/dt) + 1; 
	m_Nnodes = grid.get_nnodes(); 
	Vec_real h(m_Nnodes-1); 
	for (int i = 0; i < m_Nnodes-1; i++)
		h(i) = m_grid.get_meshsize(i);
	m_h = h; 
	
	m_Coef = Coefficient; 
	m_IncWaveFunc = IncWaveFunc; 
	
	Vec_real x(m_Nnodes); x = 0.0;  
	Vec_real x2(m_Nnodes-1); x2 = 0.0;  
	m_MassMatrix_l = x2;
	m_MassMatrix_m = x;
	m_MassMatrix_u = x2; 
	m_StiffMatrix_l = x2;
	m_StiffMatrix_m = x;
	m_StiffMatrix_u = x2;
	m_Solution_previous = x;
	m_Solution_current = x;
	m_Solution_next = x;
	m_Solution.newsize(m_NrTimeSteps, m_Nnodes); m_Solution = 0.0; 
	
	
}

WaveEq1D::~WaveEq1D() //destructor
{	m_grid = Grid1D();
	m_MaxTime = 0.0;
	m_TimeIncWaveExcited = 0.0;
	m_dt = 0.0;
	m_NrTimeSteps = 0;
	m_Nnodes = 0;

	Vec_real x(1); x = 0.0;  
	m_h = x; 
	m_Coef = x; 
	m_IncWaveFunc = x; 
	m_MassMatrix_l = x;
	m_MassMatrix_m = x;
	m_MassMatrix_u = x; 
	m_StiffMatrix_l = x;
	m_StiffMatrix_m = x;
	m_StiffMatrix_u = x;
	m_Solution_previous = x;
	m_Solution_current = x;
	m_Solution_next = x;
	m_Solution.newsize(1,1); m_Solution = 0.0; 
	
}

void WaveEq1D::MassMatrixAssembly() // compute three diagonals of the mass matrix
{	double v; 
	for (int i = 0; i < m_Nnodes-1; i++)
	{	v = (m_Coef(i) + m_Coef(i+1))*m_h(i)/12; 
		m_MassMatrix_l(i) = v;
		m_MassMatrix_u(i) = v;
	}
	m_MassMatrix_m(0) = (m_Coef(0) + m_Coef(1))*m_h(0)/6;
	for (int i = 1; i < m_Nnodes-1; i++)
	{	m_MassMatrix_m(i) = (m_Coef(i-1) + m_Coef(i))*m_h(i-1)/6 
				  + (m_Coef(i) + m_Coef(i+1))*m_h(i)/6;
	}
	m_MassMatrix_m(m_Nnodes-1) = (m_Coef(m_Nnodes-2) + m_Coef(m_Nnodes-1))*m_h(m_Nnodes-2)/6;

}

void WaveEq1D::StiffnessMatrixAssembly() // compute the three diagonals of the stiff matrix
{	for (int i = 0; i < m_Nnodes-1; i++)
	{	m_StiffMatrix_l(i) = -1/m_h(i);
		m_StiffMatrix_u(i) = -1/m_h(i);
	}
	m_StiffMatrix_m(0) = 1/m_h(0);
	for (int i = 1; i < m_Nnodes-1; i++)
	{	m_StiffMatrix_m(i) = 1/m_h(i-1) + 1/m_h(i);
	}
	m_StiffMatrix_m(m_Nnodes-1) = 1/m_h(m_Nnodes-2);	

}

void WaveEq1D::Swap() // swap the iteration: Solution_previous = Solution_current, Solution_current = Solution_next. 
{	
	for (int i = 0; i < m_Nnodes; i++)
	{	m_Solution_previous(i) = m_Solution_current(i);
		m_Solution_current(i) = m_Solution_next(i);
	}

}

void WaveEq1D::Solve_implicitFEM() //solve the equation at time;
{
	m_Solution_previous = 0.0;
	m_Solution_current = 0.0;
	
	MassMatrixAssembly(); // assemble the mass matrix
	StiffnessMatrixAssembly(); // assemble the stiffness matrix

	Vec_real A(m_Nnodes-1), B(m_Nnodes), C(m_Nnodes-1), D(m_Nnodes); // D is the right hand side vector
	A = 0.0; B = 0.0; C = 0.0; D = 0.0;
	
	//int NTimeStepsExcited = round(m_TimeIncWaveExcited/m_dt) + 1;
	int NTimeStepsExcited = m_NrTimeSteps; // TEST: use the Neumann boundary condition in the whole time 
	int n;
	A = m_MassMatrix_l/m_dt/m_dt + m_StiffMatrix_l/2; 
	B = m_MassMatrix_m/m_dt/m_dt + m_StiffMatrix_m/2; B(0) = B(0) + 1/m_dt; 
	C = m_MassMatrix_u/m_dt/m_dt + m_StiffMatrix_u/2; 
		
	for (n = 2; n < NTimeStepsExcited; n++) //solution when the incident wave is still excited
	{ 	
		D = tridiagonalmatrix_mult(m_MassMatrix_l, m_MassMatrix_m, m_MassMatrix_u, 2*m_Solution_current-m_Solution_previous)/m_dt/m_dt
		  - tridiagonalmatrix_mult(m_StiffMatrix_l, m_StiffMatrix_m, m_StiffMatrix_u, m_Solution_current)/2;
		D(m_Nnodes-1) = D(m_Nnodes-1) + m_IncWaveFunc(n);
		D(0) = D(0) + m_Solution_current(0)/m_dt; 
			
		m_Solution_next = tridiagonal_system(A,B,C,D);		
		replace_row(m_Solution,m_Solution_next,n);
		Swap();
	}
	for (n = NTimeStepsExcited; n < m_NrTimeSteps; n++) // solution after the incident stopped
	{	B(m_Nnodes-1) = B(m_Nnodes-1) + 1/m_dt; 
		
		D = tridiagonalmatrix_mult(m_MassMatrix_l, m_MassMatrix_m, m_MassMatrix_u, 2*m_Solution_current-m_Solution_previous)/m_dt/m_dt
		  - tridiagonalmatrix_mult(m_StiffMatrix_l, m_StiffMatrix_m, m_StiffMatrix_u, m_Solution_current)/2;
		D(0) = D(0) + m_Solution_current(0)/m_dt; 
		D(m_Nnodes-1) = D(m_Nnodes-1) + m_Solution_current(m_Nnodes-1)/m_dt; 

		m_Solution_next = tridiagonal_system(A,B,C,D);		
		replace_row(m_Solution,m_Solution_next,n);
		Swap();	
	}				
}


void WaveEq1D::Solve_implicitFDM_Dirichlet() //solve the equation at time;
{

	// calculate the coefficient matrix for the finite difference scheme with Dirichlet b.c. 
	// Note that we need only m_Nnodes-1 elements for A_m. 
	// However, for using with the Neumann b.c. later, we add one more elements to these vectors. 
	// These elements are not used in the Dirichlet b.c. 

	Vec_real A_l(m_Nnodes-1), A_m(m_Nnodes), A_u(m_Nnodes-1); 
	A_l = 0.0; A_m = 0.0; A_u = 0.0;

	double dt2 = m_dt*m_dt;
	int i = 0;	
	A_m(i) = dt2/m_Coef(i)/m_h(i)/m_h(i);  A_u(i) = -A_m(i); 
	for (i = 1; i < m_Nnodes-1; i++)
	{	A_l(i-1) = -2*dt2/(m_h(i-1) + m_h(i))/m_h(i-1)/m_Coef(i);
		A_u(i) = -2*dt2/(m_h(i-1) + m_h(i))/m_h(i)/m_Coef(i); // A_u(m_Nnodes-1) is not used in the Dirichlet b.c., but used in the Neumann b.c. later. 
		A_m(i) = 2*dt2/m_h(i-1)/m_h(i)/m_Coef(i);
	}
	
	// three diagonals of the coefficient matrix: 
	Vec_real A(m_Nnodes-1), B(m_Nnodes), C(m_Nnodes-1); // D is the right hand side vector

	int NTimeStepsExcited = round(m_TimeIncWaveExcited/m_dt) + 1;
	//int NTimeStepsExcited = m_NrTimeSteps; // TEST: use the Neumann boundary condition in the whole time 
	int n;

	A = A_l/2; 
	B = 1 + A_m/2; B(0) = B(0) + m_dt/m_Coef(0)/m_h(0); 
	C = A_u/2;
 
	Vec_real Up(m_Nnodes-1), Uc(m_Nnodes-1); Vec_real Un(m_Nnodes-1);
	Up = 0.0; Uc = 0.0; 
	
	m_Solution_previous = 0.0; 
	m_Solution_current = 0.0;
	
	// Dirichet boundary condition for the incident plane wave:
	Vec_real RHS(m_Nnodes-1);
	for (n = 1; n < NTimeStepsExcited; n++) //solution when the incident wave is still excited
	{ 	
		RHS = 2*Uc - Up - tridiagonalmatrix_mult(A_l, A_m, A_u, Uc)/2;
		RHS(m_Nnodes-2) = RHS(m_Nnodes-2) + (m_IncWaveFunc(n)+m_IncWaveFunc(n-1))/2*dt2/m_Coef(m_Nnodes-1)/m_h(m_Nnodes-2)/m_h(m_Nnodes-2);
		RHS(0) = RHS(0) + m_dt/m_Coef(0)/m_h(0)*m_Solution_current(0); 
		
		Un = tridiagonal_system(A,B,C,RHS);					
		m_Solution_next = vector_concate(Un,m_IncWaveFunc(n));
		replace_row(m_Solution,m_Solution_next,n);

		Up = Uc; Uc = Un; 
		Swap();
	}

	// ABC: 
	i = m_Nnodes-1;	
	A_l(i) = -dt2/m_h(i-1)/m_h(i-1)/m_Coef(i);  
	A_m(i) = -A_l(i);
	A = A_l/2; 
	B = 1 + A_m/2; 
	B(0) = B(0) + m_dt/m_Coef(0)/m_h(0); 
	B(m_Nnodes-1) = B(m_Nnodes-1) + m_dt/m_Coef(m_Nnodes-1)/m_h(m_Nnodes-2); 
	C = A_u/2;

	Vec_real RHS2(m_Nnodes); //the right hand side vector

	for (n = NTimeStepsExcited; n < m_NrTimeSteps; n++) // solution after the incident stopped
	{	
		RHS2 = (2*m_Solution_current-m_Solution_previous)
		  - tridiagonalmatrix_mult(A_l, A_m, A_u, m_Solution_current)/2;
		RHS2(0) = RHS2(0) + m_dt/m_Coef(0)/m_h(0)*m_Solution_current(0); 
		RHS2(m_Nnodes-1) = RHS2(m_Nnodes-1) + m_dt/m_Coef(m_Nnodes-1)/m_h(m_Nnodes-2)*m_Solution_current(m_Nnodes-1); 

		m_Solution_next = tridiagonal_system(A,B,C,RHS2);		
		replace_row(m_Solution,m_Solution_next,n);
		Swap();	
	}				
}


Mat_real WaveEq1D::get_solution() //load the solution to a matrix outside of the class
{
	return m_Solution;
}
















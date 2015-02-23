
#pragma once
#include <mkl.h>
#include "sys.h"
#include "FE_Storage.h"
#include "post_proc.h"
/* класс содержит все механизмы для составления СЛАУ и Решения */

#define SOL_WAIT 1
#define SOL_SOLVING 2

#define SAE_CHOLESKY 1
#define SAE_BUNCH 2
#define SAE_DSS 3
//class FE_Storage;

class Solution
{
public:
	Solution () {
		storage = NULL;
		status = SOL_WAIT;
		setInit();
		curCriteria = 0.0f;
		curLoadstep = 0;
		curIterat = 0;
		solver_type = SAE_CHOLESKY;
		if_solver_first_time = true;
	};
	void attach (FE_Storage_Interface *st);
	void run (bool multi_thread = false); // запустить решение задачи
	void stop (); // остановить решение задачи

	uint16 getqIterat ();
	uint16 getqLoadstep ();
	double getCriteria ();
	uint16 getStatus ();
	bool setqIterat (uint16 iter);
	bool setqLoadstep (uint16 ls);
	bool setCriteria (double cr);
	bool setInit ();
	
private:

	void *pt[64]; //Internal solver memory pointer pt
	bool if_solver_first_time;
	uint16 qLoadstep;	// кол-во шагов нагружений
	uint16 qIterat;		// кол-во итераций
	double Criteria;	// кинематический критерий сходимости
	uint16 status;
	uint16 solver_type;
	bool stopit;

	double curCriteria; // текущее значение критерия на данном шаге решения
	uint16 curLoadstep;
	uint16 curIterat;

	void bound (double dF_par); //применение узловых граничных условий
	uint16 solve_wrap ();
	uint16 solveSAE_Bunch ();
	uint16 solveSAE_Cholesky();
	uint16 solveSAE_DSS ();
	static void fork (void *ptr);
	FE_Storage_Interface* storage;
	void main_process (); //главный процесс решения, вызывается после отделения процесса
};

inline uint16 Solution::getqIterat ()
{
	return qIterat;
}
inline uint16 Solution::getqLoadstep ()
{
	return qLoadstep;
}
inline double Solution::getCriteria ()
{
	return Criteria;
}
inline uint16 Solution::getStatus ()
{
	return status;
}
#include "elist_proc.h"

EList_proc::EList_proc(FE_Storage *st, string str_code, uint16 gp) : Post_proc(st)
{
	code = S_UNDEF;
	calc_scheme = gp;
	name ="Elist_proc";
	char buf[200];
	if (str_code == "EX")
		code = E_X;
	else if (str_code == "EY")
		code = E_Y;
	else if (str_code == "EZ")
		code = E_Z;
	else if (str_code == "EXY")
		code = E_XY;
	else if (str_code == "EXZ")
		code = E_XZ;
	else if (str_code == "EYZ")
		code = E_YZ;
	else if (str_code == "ER")
		code = E_R;
	else if (str_code == "ET")
		code = E_T;
	else if (str_code == "ERZ")
		code = E_RZ;
	else if (str_code == "ERT")
		code = E_RT;
	else if (str_code == "EZT")
		code = E_ZT;
	else if (str_code == "EVOL")
		code = E_VOL;
	else if (str_code == "SX")
		code = S_X;
	else if (str_code == "SY")
		code = S_Y;
	else if (str_code == "SZ")
		code = S_Z;
	else if (str_code == "SXY")
		code = S_XY;
	else if (str_code == "SXZ")
		code = S_XZ;
	else if (str_code == "SYZ")
		code = S_YZ;
	else if (str_code == "SR")
		code = S_R;
	else if (str_code == "ST")
		code = S_T;
	else if (str_code == "SRZ")
		code = S_RZ;
	else if (str_code == "SRT")
		code = S_RT;
	else if (str_code == "SZT")
		code = S_ZT;
	else if (str_code == "SP")
		code = S_P;
	else {
		warning ("Unkown label %s", str_code.c_str());
		failed = true;
		return;
	}
	name += "("+str_code+")";
	sprintf_s(buf, 200, "POST%d %s",nPost_proc, name.c_str());
	file_name = buf;
	elems.clear();
	every_n_ls = 1;
	active = true;
}
EList_proc::~EList_proc()
{

}
void EList_proc::pre (uint16 qLoadstep)
{
	ofstream file(file_name, ios::trunc);
	file.close();
}
void EList_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	if (!active) return;
	//проверка: подходящая ли конфигурация
	if (calc_scheme != GP_MEAN && calc_scheme >= _4N_2D_Form::nInt*_4N_2D_Form::nInt)
	{
		warning("EList_proc::process: calc_scheme isn't appropriate for the current integration scheme");
		return;
	}
	//if (остаток от деления curLoadstep на every_n_ls)
	ofstream file(file_name, ios::app);
	file << "LoadStep " << curLoadstep << endl;
	double val;
	uint16 g_start;
	uint16 g_end;
	if (calc_scheme == GP_MEAN)
	{
		g_start = 0;
		g_end = _4N_2D_Form::nInt*_4N_2D_Form::nInt;
	} else
	{
		g_start = calc_scheme;
		g_end = g_start+1;
	}
	for (uint32 el=0; el < elems.size(); el++)
	{
		val=0.0f;
		for (uint16 j = g_start; j < g_end ; j++)
		{
			val += storage->el_form[elems[el]]->getComponent(j,code, elems[el], storage);
		}
		val = val/(float)(g_end-g_start);
		file << val << ";";
	}
	file << endl << endl;
	file.close();
}
void EList_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{

}

string EList_proc::getStatus()
{
	return Post_proc::getStatus();
}

#include "lambda_proc.h"

lambda_proc::lambda_proc(FE_Storage *st, string _fileName) : Post_proc(st) {
	name = "lambda_processor";
	file_name = _fileName;
}

lambda_proc::~lambda_proc() {
}

void lambda_proc::pre(uint16 qLoadstep) {
}

void lambda_proc::process (uint16 curLoadstep, uint16 qLoadstep)
{
	string cur_fn = file_name +  IntToStr(curLoadstep) + ".txt";
	ofstream file(cur_fn, ios::trunc);
  file << "element\tvolume\tlambda1\tlambda2\tlambda3" << endl;
	for (uint32 i=1; i <= storage->getNumElement(); i++) {
    MIXED_8N_3D_P0* el = dynamic_cast<MIXED_8N_3D_P0*> (&storage->getElement(i));
    if (el == NULL) {
      error("lambda_proc::process: Can't cast to el");
    }
    double dWtSym = 0.0;
    double dWt = 0.0;
    double lambda1 = 0.0;
    double lambda2 = 0.0;
    double lambda3 = 0.0;
    for (uint16 nPoint = 0; nPoint < npow(Element::n_int(),Element::n_dim()); nPoint ++) {
      dWt = el->g_weight(nPoint);
      dWtSym += dWt;

      lambda1 += sqrt(el->getComponent(nPoint,E_1,i,storage))*dWt;
      lambda2 += sqrt(el->getComponent(nPoint,E_2,i,storage))*dWt;
      lambda3 += sqrt(el->getComponent(nPoint,E_3,i,storage))*dWt;
    }
    lambda1 = lambda1 / dWtSym;
    lambda2 = lambda2 / dWtSym;
    lambda3 = lambda3 / dWtSym;
    file << i << "\t" << dWtSym << "\t" << lambda1 << "\t" << lambda2 << "\t" << lambda3 << endl;
  }
	file.close();
}

void lambda_proc::post (uint16 curLoadstep, uint16 qLoadstep)
{
}


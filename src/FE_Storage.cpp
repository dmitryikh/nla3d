
#include "FE_Storage.h"


//bool FE_Storage::save_mdl (const char *filename)
//{
//	ofstream file(filename);
//	if (!file)
//	{
//		warning("FE_Storage::save_mdl: Can't create save fale `%s`",filename);
//		return false;
//	}
//	file << "2D" << " ";
//	if (fem_type==FEM_PLANE)
//		file << "PLANE ";
//	else if (fem_type==FEM_AXYM)
//		file << "AXYM ";
//	else
//	{
//		warning("FE_Storage::save_mdl: unknown type of fem_type = %d",fem_type);
//		file.close();
//		return false;
//	}
//	file << nn << " " << en << " " << bounds.size() << endl;
//	for (uint32 i=0; i < nn; i++)
//		file << nodes[i].toString() << endl;
//	for (uint32 i=0; i < en; i++)
//		file << elements[i].toString() << endl;
//	for (uint32 i=0; i < bounds.size(); i++)
//		file << bounds[i].toString() << endl;
//	if (material)
//		file << material->toString() << endl;
//	file.close();
//	return true;
//}

//bool FE_Storage::load_mdl (const char *filename)
//{
//	//удаляем всё
//	//TODO: залочить доступ к storage
//	if (material) delete material;
//	elements.clear();
//	nodes.clear();
//	bounds.clear();
//	nn=0;
//	en=0;
//	nDOFs=0;
//
//	ifstream file(filename);
//	if (!file)
//	{
//		warning("FE_Storage::load_mdl: Can't open file `%s`",filename);
//		return false;
//	}
//
//	string str;
//	file >> str;
//	if (str != "2D")
//	{
//		warning("FE_Storage::load_mdl: unknown dimension of solution (%s)", str.c_str());
//		file.close();
//		return false;
//	}
//	file >> str;
//	if (str=="PLANE")
//		fem_type = FEM_PLANE;
//	else if (str=="AXYM")
//		fem_type = FEM_AXYM;
//	else
//	{
//		warning("FE_Storage::load_mdl: unknown type of analysis (%s)",str);
//		file.close();
//		return false;
//	}
//	file >> nn;
//	if (nn < 1)
//	{
//		warning("FE_Storage::load_mdl: nodes number is invalid (%d)",nn);
//		file.close();
//		return false;
//	}
//	file >> en;
//	if (en <1)
//	{
//		warning("FE_Storage::load_mdl: elements number is invalid (%d)",en);
//		file.close();
//		return false;
//	}
//	uint32 bcn;
//	file >> bcn;
//	if (bcn < 1)
//	{
//		warning("FE_Storage::load_mdl: boundary conditions number is invalid (%d)",bcn);
//		file.close();
//		return false;
//	}
//	nodes.assign(nn, Node());
//	elements.assign(en, Element4_2D());	
//	bounds.assign(bcn, Bound());
//	for (uint32 i=0; i < nn; i++)
//		nodes[i].read_from_stream(file);
//	for (uint32 i=0; i < en; i++)
//		elements[i].read_from_stream(file);
//	for (uint32 i=0; i < bcn; i++)
//		bounds[i].read_from_stream(file);
//	file >> str;
//	if (str == "Material_Comp_Neo_Hookean")
//	{
//		material = new Material_Comp_Neo_Hookean();
//		file >> material->E;
//		file >> material->mu;
//	}
//	else if (str == "Material_Hookean")
//	{
//			material = new Material_Hookean();
//			file >> material->E;
//			file >> material->mu;
//	}
//	status = ST_LOADED;
//	nDOFs = nn*Node::nDOF;
//	file.close();
//	return true;
//}
//


//void FE_Storage::makestripmodel(double B, double H, uint16 nB, uint16 nH, double alpha, double r, double delta)
//{
//	uint32 nN = (nB+1)*(nH+1);
//	uint32 nE = nB*nH;
//	nodes.clear();
//	nodes.assign(nN, Node());
//	elements.clear();
//	elements.assign(nE, Element4_2D());
//	bounds.clear();
//	en = nE;
//	nn = nN;
//	nDOFs=nn*Node::nDOF;
//	Force_proc *fproc = new Force_proc(this);
//	EList_proc *el_proc1 = new EList_proc(this, "EVOL", 4);
//	EList_proc *el_proc2 = new EList_proc(this, "EY", 4);
//	EList_proc *el_proc3 = new EList_proc(this, "EXY", 4);
//	EList_proc *el_proc4 = new EList_proc(this, "SX", GP_MEAN);
//	EList_proc *el_proc5 = new EList_proc(this, "SY", GP_MEAN);
//	EList_proc *el_proc6 = new EList_proc(this, "SXY", GP_MEAN);
//	for (uint16 i=0; i < nH+1; i++)
//	{
//		for (uint16 j=0; j < nB+1; j++)
//		{
//			nodes[i*(nB+1)+j].coord[0] = r+(H/nH)*i*tan(M_PI/2-alpha)+j*B/nB;
//			nodes[i*(nB+1)+j].coord[1] = (H/nH)*i;
//			if (i==nH)
//			{
//				fproc->nodes.push_back(i*(nB+1)+j+1);
//				fproc->dofs.push_back(1);
//			}
//			if (i>0 && j>0)
//			{
//				elements[(i-1)*nB+j-1].nodes[0]=(i-1)*(nB+1)+j-1+1;
//				elements[(i-1)*nB+j-1].nodes[3]=(i)*(nB+1)+j-1+1;
//				elements[(i-1)*nB+j-1].nodes[2]=(i)*(nB+1)+j+1;
//				elements[(i-1)*nB+j-1].nodes[1]=(i-1)*(nB+1)+j+1;
//				if ((int)((nH)/4) == i)
//				{
//					el_proc1->elems.push_back((i-1)*nB+j-1);
//					el_proc2->elems.push_back((i-1)*nB+j-1);
//					el_proc3->elems.push_back((i-1)*nB+j-1);
//				}
//				if ((int)((nH)/2) == i)
//				{
//					el_proc4->elems.push_back((i-1)*nB+j-1);
//					el_proc5->elems.push_back((i-1)*nB+j-1);
//					el_proc6->elems.push_back((i-1)*nB+j-1);
//				}
//			}
//		}
//	}
//	// boundary constrains
//	for (uint16 j=0; j<nB+1; j++)
//	{
//		Bound bd;
//		bd.key = D_UX;
//		bd.value = 0.0f;
//		bd.node = j+1;
//		bounds.push_back(bd);
//		bd.key = D_UY;
//		bounds.push_back(bd);
//		bd.key = D_UX;
//		bd.node = ((nH)*(nB+1))+j+1;
//		bounds.push_back(bd);
//		bd.key = D_UY;
//		bd.value = -H*delta;
//		bounds.push_back(bd);
//	}
//	status=ST_LOADED;
//}
//
//void FE_Storage::makecircmodel (double r1, double r2, double fi1, double fi2, uint16 nR, uint16 nFi, double P)
//{
//	uint32 nN = (nR+1)*(nFi+1);
//	uint32 nE = nR*nFi;
//	nodes.clear();
//	nodes.assign(nN, Node());
//	elements.clear();
//	elements.assign(nE, Element4_2D());
//	bounds.clear();
//	en = nE;
//	nn = nN;
//	nDOFs=nn*Node::nDOF;
//	double dFi = (fi2 - fi1)/(nFi);
//	double dR = (r2-r1)/(nR);
//	double dP = P/(nR+1);
//	EList_proc *el_proc1 = new EList_proc(this, "SX", GP_MEAN);
//	EList_proc *el_proc2 = new EList_proc(this, "SY", GP_MEAN);
//	EList_proc *el_proc3 = new EList_proc(this, "SXY", GP_MEAN);
//	NList_proc *n_proc4 = new NList_proc(this, "UX");
//	for (uint16 i=0; i < nFi+1; i++)
//	{
//		for (uint16 j=0; j < nR+1; j++)
//		{
//			nodes[i*(nR+1)+j].coord[0] = (r1+dR*j)*cos(fi1+dFi*i);
//			nodes[i*(nR+1)+j].coord[1] = (r1+dR*j)*sin(fi1+dFi*i);
//			if (i>0 && j>0)
//			{
//				elements[(i-1)*nR+j-1].nodes[0]=(i-1)*(nR+1)+j-1+1;
//				elements[(i-1)*nR+j-1].nodes[3]=(i)*(nR+1)+j-1+1; //
//				elements[(i-1)*nR+j-1].nodes[2]=(i)*(nR+1)+j+1; //
//				elements[(i-1)*nR+j-1].nodes[1]=(i-1)*(nR+1)+j+1;
//				if ((int)((nFi)/3) == i)
//				{
//					el_proc1->elems.push_back((i-1)*nR+j-1);
//					el_proc2->elems.push_back((i-1)*nR+j-1);
//					el_proc3->elems.push_back((i-1)*nR+j-1);
//				}
//			}
//			if ((int)((nFi)/3) == i)
//			{
//				n_proc4->nodes.push_back(i*(nR+1)+j);
//			}
//		}
//	}
//	// boundary constrains
//	uint16 i_zad = nFi;
//	uint16 i_p = 0;
//	for (uint16 j=0; j<nR+1; j++)
//	{
//		Bound bd;
//		// заделка
//		bd.key = D_UX;
//		bd.value = 0.0f;
//		bd.node = i_zad*(nR+1)+j+1;
//		bounds.push_back(bd);
//		bd.key = D_UY;
//		bounds.push_back(bd);
//		// сила
//		bd.key = F_X;
//		bd.node = i_p*(nR+1)+j+1;
//		bd.value = -dP*cos(fi1);
//		bounds.push_back(bd);
//		bd.key = F_Y;
//		bd.value = -dP*sin(fi1);
//		bounds.push_back(bd);
//	}
//	status=ST_LOADED;
//}


//
//double FE_Storage::node_q(uint32 n, uint16 dof)
//{
//	assert(dof >= 0 && dof < Node::nDOF);
//	//dof: x - 0, y - 1, z - 2
//	// n: от 1 
//	return (*vecQsum)[(n-1)*Node::nDOF+dof+1];
//}

//void FE_Storage::save_to_fortran (const char *filename)
//{
//	ofstream file(filename);
//	file << en << " " << nn << " " << 2 << " " << 4 << endl;
//	for (uint32 i=0; i < nn; i++)
//	{
//		file << i+1 << " " << nodes[i].coord[0] << " " << nodes[i].coord[1] << endl;
//	}
//	for (uint32 i=0; i < en; i++)
//	{
//		file << i+1 << " " << elements[i].nodes[0] << " " << elements[i].nodes[1] << " " << elements[i].nodes[2] << " " << elements[i].nodes[3] << endl;
//	}
//	file << bounds.size() << endl;
//	for (uint32 i=0; i < bounds.size(); i++)
//	{
//		file << bounds[i].node << " ";
//		if (bounds[i].key == D_UX) file << "10 " << bounds[i].value << " 0.0" << endl;
//		else file << "01 0.0 " <<bounds[i].value << endl;
//	}
//	file << "6.0  0.45" << endl << "3   8" << endl << 0 << endl << 6;
//	file.close();
//}


// read Ansys Mechanical APDL *.cdb file. Nodes, Elements, Displacement BC and MPC (Constraint equations) is supported
// read_ans_data repcales storage's mesh.
bool read_ans_data(const char *filename, FE_Storage_Interface *storage)
{
	uint32 n_number, en;
	ifstream file(filename);
	if (!file)
	{
		warning("read_ans_data: Can't open input file `%s`",filename);
		return false;
	}
	storage->clearMesh();
	char buf[1024]="";
	while (file.getline(buf, 1024))
	{
		vector<string> vec = read_tokens(buf);
		if (vec[0].compare("NBLOCK") == 0)
		{
    //NBLOCK,6,SOLID,     9355,     9355
    //(3i9,6e20.13)
    //        1        0        0 7.0785325971794E+01 6.5691449317818E+01-3.6714639015390E+01
			uint32 max_n_number = atoi(vec[6].c_str());
			n_number= atoi(vec[6].c_str());
      if (max_n_number != n_number) {
        warning("read_ans_data: NBLOCK: maximum node number is %d, but number of nodes is %s. Note that nla3d needs compressed numbering for nodes and elements", max_n_number, n_number );
      }
			storage->nodes_reassign(n_number);
			file.getline(buf, 1024);
      string buf_str(buf);
      // we need to take a format of columns "3i9"
      uint16 start = buf_str.find("i")+1;
      uint16 stop = buf_str.find(",");
      string tmp;
      tmp.assign((char*) (buf + start), stop-start);
			uint16 frmt = atoi(tmp.c_str());

      //TODO: as we need numbering from 1/0 to n_number, here we can check that number of nodes and theirs id are in a row
			for (uint32 i=1; i<=n_number; i++)
			{
				file.getline(buf, 1024);
				uint16 len=strlen(buf);
				for (uint16 j=0; j<3;j++)
					if (len>=3*frmt+20*(j+1))
            //note that last column in NBLOCK table could be avoided if Z=0.0
            //but storage->nodes_reassign(n_number) initialize the node table with (0,0,0)
						storage->getNode(i).pos[j] = atof(string((char*) (buf+3*frmt+20*j),20).c_str());
			}
		}//NBLOCK
		else if (vec[0].find("EBLOCK")!=vec[0].npos)
		{  
      //EBLOCK,19,SOLID,      7024,      7024
      //(19i9)
			en=atoi(vec[6].c_str());
			storage->elements_reassign(en);
			file.getline(buf, 1024);
      // we need to take a format of columns "3i9"
      // in Ansys 12 here is 8 symbols per number (19i8), but in ansys 15 (19i9) is used. 
      string buf_str(buf);
      uint16 start = buf_str.find("i")+1;
      uint16 stop = buf_str.find(")");
      string tmp;
      tmp.assign((char*) (buf + start), stop-start);
			uint16 frmt = atoi(tmp.c_str()); 
			for (uint32 i=1; i<=en; i++)
			{
				file.getline(buf, 1024);
				uint16 len=strlen(buf);
        //TODO: check that n_nodes and provided number of nodes are the same
        if (len != 11*frmt+frmt*Element::n_nodes())
          warning("read_ans_data: in EBLOCK for element %d the number of nodes provided is not equal to %d", i, Element::n_nodes());
				for (uint16 j=0; j<Element::n_nodes();j++)
					if (len>=11*frmt+frmt*(j+1))
						storage->getElement(i).node_num(j) = atoi(string((char*) (buf+11*frmt+frmt*j),frmt).c_str());
			}
		}//EBLOCK
		else if (vec[0].find('D')!=vec[0].npos && vec[0].length()==1)
		{
				BC_dof_constraint bnd;
				bnd.node = atoi(vec[2].c_str());
				bnd.value = atof(vec[6].c_str());
        bnd.node_dof = str2dof(vec[4]);
        storage->add_bounds(bnd);
		}//D
    else if (vec[0].compare("CE") == 0)
    {
      //How MPC looks like this in cdb file:
      //CE,R5.0,DEFI,       2,       1,  0.00000000    
      //CE,R5.0,NODE,      1700,UX  ,  1.00000000    ,      1700,UZ  ,  1.00000000  
      BC_MPC mpc;
      mpc.b = atoi(vec[10].c_str()); //rhs of MPC equation
      uint16 n_terms = atoi(vec[6].c_str()); //number of terms in equation
      //debug("MPC link: %d terms, b = %f", n_terms, mpc.b);
      while (n_terms > 0)
      {
        file.getline(buf, 1024);
        vector<string> vec = read_tokens(buf);
        uint16 place = 6;
        for (int i=0; i < max(n_terms, 2); i++) 
        {
          uint32 node = atoi(vec[place].c_str());
          uint16 dof = str2dof(vec[place+2]);
          double coef = atof(vec[place+4].c_str());
          //debug("%d term: node %d, dof %d, coef = %f", i, node, dof, coef);
          mpc.eq.push_back(MPC_token(node,dof,coef));
          place += 6;
          n_terms--;
        }
      }
			storage->add_bounds(mpc);
    }//CE (MPC)
    //TODO: add FX FY FZ
	}
	file.close();
	storage->setStatus(ST_LOADED);
	return true;
}

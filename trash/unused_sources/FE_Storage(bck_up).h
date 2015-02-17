#pragma once
#include <string>
#include "sys.h"
#include "Material.h"
#include "Element.h"
#include "math\Mat_Band.h"
#include "post_proc.h"
 
//pre-defines
class Element;
class Node;
class MIXED_8N_3D_Form;
class Post_proc;

// class Bound represents boundary condition on the node (disp & force BC)
// DIM INDEPENDENT CALSS
// TODO: add force support
struct Bound 
{
	uint32 node;
	uint16 key;
	double value;
	void display(uint32 bn)
	{
		string k="UNKNOWN";
		switch (key)
		{
			case D_UX: k = "UX";
				break;
			case D_UY: k = "UY";
				break;
			case D_UZ: k = "UZ";
				break;
		}
		echo("B %d: %d\t%s\t%f",bn, node, k.c_str(),value);
	}
	string toString()
	{
		string str;
		char buff[100];
		switch (key)
		{
			case D_UX: str="UX";
				break;
			case D_UY: str="UY";
				 break;
			case D_UZ: str="UZ";
				 break;
			default:
				warning("Bound::toString: unknown bound key %d", key);
		}
		sprintf_s(buff,100," %d %f",node, value);
		str+=buff;
		return str;
	}
	void read_from_stream (istream &file)
	{
		string str;
		file >> str;
		if (str=="UX")
			key = D_UX;
		else if (str=="UY")
			key = D_UY;
		else if (str=="UZ")
			key = D_UZ;
		else
		{
			//TODO: error
			warning("Bound::read_from_stream: unknown bound key %d", key);
		}
		file >> node;
		file >> value;
	}
};

#define ST_INIT 1
#define ST_LOADED 2
#define ST_SOLVED 3

// abstract class to provide interface to Storage
class FE_Storage_Interface
{
public:
	virtual double& getKij(uint32 nodei, uint16 dofi, uint32 nodej, uint16 dofj) = 0;
	virtual double& getFi(uint32 nodei, uint16 dofi) = 0;
	virtual void clearK() = 0;
	virtual void clearF() = 0;
	virtual	Mat_Band_cm& getMatK() = 0;
	virtual Vec_Long& getVecF () = 0;
	virtual Vec_Long& getVecdQ() = 0;
	virtual Vec_Long& getVecQsum() = 0;

	virtual Material& getMaterial() = 0;
	virtual Node& getNode(uint32 _nn) = 0;
	virtual Element& getElement(uint32 _en) = 0;
	virtual Post_proc& getPostProc(uint32 _np) = 0;
	virtual Bound& getBound (uint32 _nb) = 0;

	virtual uint32 get_nBand () = 0;
	virtual uint32 getNumDofs () = 0;
	virtual uint32 getNumNode () = 0;
	virtual uint32 getNumElement () = 0;
	virtual uint32 getNumPostProc () = 0;
	virtual uint32 getNumBound() = 0;
	virtual uint16 getStatus() = 0;
	virtual void setStatus(uint16 st) = 0;

	virtual void add_Bounds (Bound &bnd) = 0;
	virtual uint16 add_post_proc (Post_proc *pp) = 0;

	virtual void clearMesh () = 0;
	virtual void nodes_reassign(uint32 nn) = 0; // dont forget 	nDOFs = nn*Node::nDOF;
	virtual void elements_reassign(uint32 en) = 0;

	virtual bool prepare_for_solution ()=0;

	virtual void element_nodes(uint32 el, Node** node_ptr) = 0;
	virtual void get_q_e(uint32 el, double* ptr) = 0;
	virtual void get_dq_e(uint32 el, double* ptr) = 0;
	virtual void get_q_n(uint32 n, double* ptr) = 0;
	virtual double get_qi_n(uint32 n, uint16 dof) = 0;
	virtual void get_node_pos(uint32 n, double* ptr, bool def = false) = 0; //def. or initial pos

};

// class FE_Storage - heart of the NLA programm, it contains all data to keep in memory: 
//				nodes, elements, solution staff (global matrixes and vectors), processors.
//TODO: add mutexes to prevent multithreads errors
template <typename el_type>
class FE_Storage : public FE_Storage_Interface
{
public:
	FE_Storage();
	~FE_Storage();
	el_type dummy; // пустышка, для разных нужд
	vector<el_type> elements;
	vector<Node> nodes;
	vector<Bound> bounds;
	vector<Post_proc*> post_procs;

	uint32 nn;
	uint32 en;
	uint32 nDOFs; //общее число степеней свободы системы
	uint32 nBand;

	Material *material;
	Mat_Band_cm *matK; //глобальная матрица жесткости конструкции (column oriented для MKL решателей)
	Vec_Long *vecF; //вектор нагрузки от накопленных деформаций
	Vec_Long *vecQsum; //вектор накопленных перемещений
	Vec_Long *vecdQ; // приращения перемещений

	//status
	uint16 status;

	//bool save_mdl (const char *filename);  // convert to 3d
	//bool load_mdl (const char *filename);  // convert to 3d
	//void add_mdl (char *filename);
	//void save_to_fortran (const char *filename);  // convert to 3d
	
	//void makestripmodel (double B, double H, uint16 nB, uint16 nH, double alpha, double r, double delta);  //convert to 3d
	//void makecircmodel (double r1, double r2, double fi1, double fi2, uint16 nR, uint16 nFi, double P);

	//реализация интерфейса FE_Storage_Interface
	double& getKij(uint32 nodei, uint16 dofi, uint32 nodej, uint16 dofj);
	double& getFi(uint32 nodei, uint16 dofi);
	void clearK();
	void clearF();
	Mat_Band_cm& getMatK();
	Vec_Long& getVecF ();
	Vec_Long& getVecdQ();
	Vec_Long& getVecQsum();
	Material& getMaterial();
	Node& getNode(uint32 _nn);
	Element& getElement(uint32 _en);
	Post_proc& getPostProc(uint32 _np);
	Bound& getBound (uint32 _nb);
	uint32 get_nBand ();
	uint32 getNumDofs () {
		return nDOFs; 
	};
	uint32 getNumNode () {
		return nn;
	};
	uint32 getNumElement () {
		return en;
	};
	uint32 getNumPostProc () {
		return post_procs.size();
	};
	uint32 getNumBound() {
		return bounds.size();
	};
	uint16 getStatus() {
		return status;
	};
	void setStatus(uint16 st) {
		status = st;
	};
	uint16 add_post_proc (Post_proc *pp);
	void add_Bounds (Bound &bnd);
	void clearMesh ();
	void nodes_reassign(uint32 nn); // dont forget 	nDOFs = nn*Node::nDOF;
	void elements_reassign(uint32 en);

	bool prepare_for_solution ();

	void element_nodes(uint32 el, Node** node_ptr);
	void get_q_e(uint32 el, double* ptr);
	void get_dq_e(uint32 el, double* ptr);
	void get_q_n(uint32 n, double* ptr);
	double get_qi_n(uint32 n, uint16 dof);
	void get_node_pos(uint32 n, double* ptr, bool def = false); //def. or initial pos
private:
	void count_nBand ();
};


template <typename el_type>
FE_Storage<el_type>::FE_Storage() : dummy()
{
	matK = NULL;
	vecF = NULL;
	vecQsum = NULL;
	vecdQ = NULL;
	material = NULL;
	nn = 0;
	en = 0;
	nDOFs = 0;
	nBand = 0;
	status = ST_INIT;
};
template <typename el_type>
FE_Storage<el_type>::~FE_Storage ()
{
	status = ST_INIT; //TODO: check
	for (uint16 i=0; i < post_procs.size(); i++)
		delete post_procs[i];
	if (matK) delete matK;
	if (vecF) delete vecF;
	if (vecQsum) delete vecQsum;
	if (vecdQ) delete vecdQ;
	if (material) delete material;
	elements.clear();
	nodes.clear();
	bounds.clear();
}

template <typename el_type>
bool FE_Storage<el_type>::prepare_for_solution ()
{
	if (status < ST_LOADED)
	{
		warning ("FE_Storage::prepare_for_solution: model isn't loaded");
		return false;
	}
	nDOFs=nn*Node::nDOF();
	get_nBand();//подсчитываем ширину полуленты
	echolog("DoFs = %d, nBand = %d", nDOFs, nBand);
	if (!material)
	{
		warning("FE_Storage::prepare_for_solution: material isn't defined");
		return false;
	}
	if (matK) 
	{
		warning("FE_Storage::prepare_for_solution: matK isn't empty yet. Deleting old.");
		delete matK;
	}
	matK = new Mat_Band_cm(nDOFs, nBand);
	
	if (vecF)
	{
		warning("FE_Storage::prepare_for_solution: vecF isn't empty yet. Deleting old.");
		delete vecF;
	}
	vecF = new Vec_Long(nDOFs);

	if (vecQsum)
	{
		warning("FE_Storage::prepare_for_solution: vecQsum isn't empty yet. Deleting old.");
		delete vecQsum;
	}
	vecQsum = new Vec_Long(nDOFs);

	if (vecdQ)
	{
		warning("FE_Storage::prepare_for_solution: vecdQ isn't empty yet. Deleting old.");
		delete vecdQ;
	}
	vecdQ = new Vec_Long(nDOFs);
	return true;
}

template <typename el_type>
uint32 FE_Storage<el_type>::get_nBand ()
{
	if (!nBand) count_nBand();
	return nBand;
}

template <typename el_type>
void FE_Storage<el_type>::count_nBand ()
{
	uint32 min, max, MAX;
	MAX = 0;
	for (uint32 i=0; i<en; i++)
	{
		
		min = elements[i].node(0);
		max = min;
		for (uint16 j=1; j < Element::nNodes(); j++) 
		{
			if (elements[i].node(j) < min)
				min = elements[i].node(j);
			if (elements[i].node(j) > max)
				max = elements[i].node(j);
		}
		if ((max-min) > MAX)
			MAX = max - min;
	}
	nBand = (MAX+1)*Node::nDOF(); // полуширина ленты с центральной диагональю
}

//функция возвращает ссылку на значение, которое хранится в глобальной матрице жосткости по адресу 
// строка: dofi-тая степень свободы nodei-того узла 
// столбец: dofj-тая степень свободы nodej-того узла 
// узлы нумеруются с индекса 1, dofi=[0;Node::nDOF()-1];
// проверка на ошибки индексов (выход за пределы) производит сама структура хранения глобальной МЖ
// TODO: CHECK
template<typename el_type>
double& FE_Storage<el_type>::getKij(uint32 nodei, uint16 dofi, uint32 nodej, uint16 dofj)
{
	assert(matK);
	uint32 row = Node::nDOF()*(nodei-1) + dofi + 1;
	uint32 col = Node::nDOF()*(nodej-1) + dofj + 1;
	return (*matK)[row][col];
}
//getFi(..) возвращает ссылку на ячейку глобального вектора сил
//по адресу dofi-тая степень свободы nodei-того узла
// узлы нумеруются с индекса 1, dofi=[0;Node::nDOF()-1];
// TODO: CHECK
template<typename el_type>
double& FE_Storage<el_type>::getFi(uint32 nodei, uint16 dofi)
{
	assert(vecF);
	uint32 row = Node::nDOF()*(nodei-1) + dofi + 1;
	return (*vecF)[row];
}
// clearK() очищает матрицу жесткости
template<typename el_type>
void FE_Storage<el_type>::clearK()
{
	assert(matK);
	matK->zero();
}
// clearF() очищает глобальный вектор сил
template<typename el_type>
void FE_Storage<el_type>::clearF()
{
	assert(vecF);
	vecF->zero();
}
// getMatK возвращает ссылку на структуру данных глобальной матрицы жесткости
template<typename el_type>
Mat_Band_cm& FE_Storage<el_type>::getMatK()
{
	assert(matK);
	return *matK;
}

//getVecF()
template<typename el_type>
Vec_Long& FE_Storage<el_type>::getVecF ()
{
	assert(vecF);
	return *vecF;
}

//getVecdQ()
template<typename el_type>
Vec_Long& FE_Storage<el_type>::getVecdQ ()
{
	assert(vecdQ);
	return *vecdQ;
}

//getVecQ()
template<typename el_type>
Vec_Long& FE_Storage<el_type>::getVecQsum ()
{
	assert(vecQsum);
	return *vecQsum;
}

//getMaterial()
template<typename el_type>
Material& FE_Storage<el_type>::getMaterial()
{
	assert(material);
	return *material;
}

//getNode(nn)
//нумерация с 1
template<typename el_type>
Node& FE_Storage<el_type>::getNode(uint32 _nn)
{
	assert(_nn <= nn);
	return nodes[_nn-1];
}

//getElement(_en)
//нумерация с 1
template<typename el_type>
Element& FE_Storage<el_type>::getElement(uint32 _en)
{
	assert(_en <= en);
	return elements[_en-1];
}

//getPostProc(_np)
//нумерация с 0
template<typename el_type>
Post_proc& FE_Storage<el_type>::getPostProc(uint32 _np)
{
	assert(_np < getNumPostProc());
	return *post_procs[_np];
}

//getBound(_nb)
//нумерация с 0
template<typename el_type>
Bound& FE_Storage<el_type>::getBound (uint32 _nb)
{
	assert(_nb < getNumBound());
	return bounds[_nb];
}

//add_post_proc(*pp)
template<typename el_type>
uint16 FE_Storage<el_type>::add_post_proc (Post_proc *pp)
{
	assert(pp);
	uint16 num = this->post_procs.size()+1;
	pp->nPost_proc = num;
	post_procs.push_back(pp);
	return num;
}

//add_Bounds
template<typename el_type>
void FE_Storage<el_type>::add_Bounds (Bound &bnd)
{
	bounds.push_back(bnd);
}

//clearMesh() функция очищает таблицу узлов, элементов, ГУ
template<typename el_type>
void FE_Storage<el_type>::clearMesh ()
{
	status = ST_INIT; //TODO продумать
	elements.clear();
	nodes.clear();
	bounds.clear();
	nn = 0;
	en = 0;
	nDOFs = 0;
	nBand = 0;
}

//nodes_reassign(_nn)
template<typename el_type>
void FE_Storage<el_type>::nodes_reassign(uint32 _nn)
{
	nodes.clear();
	nn = _nn;
	nDOFs = nn*Node::nDOF();
	nodes.assign(_nn, Node()); // TODO: тут бы catch на возможность выделения памяти
}

//elements_reassign(_en)
template<typename el_type>
void FE_Storage<el_type>::elements_reassign(uint32 _en)
{
	elements.clear();
	en = _en;
	elements.assign(_en, el_type());
	nBand = 0;
}

// element_nodes(el, node_ptr), вызывающая сторона должна предоставить массив 
// указателей Node* на >= Element::nNodes() элементов
// el начинается с 1
template<typename el_type>
void FE_Storage<el_type>::element_nodes(uint32 el, Node** node_ptr)
{
	assert(el <= en);
	for (uint16 i=0; i<Element::nNodes(); i++)
		node_ptr[i] = & nodes[elements[el-1].node(i)-1]; //TODO: CHECK
}

// get_q_e(el, ptr) функция возвращает вектор узловых степеней свободы элемента,
// вызывающая сторона должна предоставить массив ptr размерностью Element::nNodes()*Node::nDOF() 
// указателей Node* на >= Element::nNodes() элементов
// el начинается с 1
template<typename el_type>
void FE_Storage<el_type>::get_q_e(uint32 el, double* ptr)
{
	assert(el <= en);
	assert(vecQsum);
	for (uint16 i=0; i<Element::nNodes(); i++)
		for (uint16 j=0; j<Node::nDOF(); j++)
			ptr[i*Node::nDOF()+j] = (*vecQsum)[Node::nDOF()*(elements[el-1].node(i)-1) + j + 1]; //TODO: CHECK
}
// get_dq_e см. выше
template<typename el_type>
void FE_Storage<el_type>::get_dq_e(uint32 el, double* ptr)
{
	assert(el <= en);
	assert(vecdQ);
	for (uint16 i=0; i<Element::nNodes(); i++)
		for (uint16 j=0; j<Node::nDOF(); j++)
			ptr[i*Node::nDOF()+j] = (*vecdQ)[Node::nDOF()*(elements[el-1].node(i)-1) + j + 1]; //TODO: CHECK
}

//get_q_n(n, ptr)
// n начинается с 1
template<typename el_type>
void FE_Storage<el_type>::get_q_n(uint32 n, double* ptr)
{
	assert(n > 0 && n < nn);
	assert(vecQsum);
	for (uint16 j=0; j<Node::nDOF(); j++)
		ptr[j] = (*vecQsum)[Node::nDOF()*(n-1)+j+1]; //TODO: CHECK
}

//get_qi_n(n, dof)
// n начинается с 1
template<typename el_type>
double FE_Storage<el_type>::get_qi_n(uint32 n, uint16 dof)
{
	assert(vecQsum);
	assert(n > 0 && n <= nn);
	assert(dof < Node::nDOF());
	return (*vecQsum)[Node::nDOF()*(n-1)+dof+1];
}

//void get_node_pos(uint32 n, double* ptr, bool def = false) double массим на 3 элемента!
//n с 1
template<typename el_type>
void FE_Storage<el_type>::get_node_pos(uint32 n, double* ptr, bool def = false)
{
	assert(n > 0 && n <= nn);
	for (uint16 i=0; i<3; i++)
		ptr[i] = nodes[n-1].vecXYZ[i];
	if (def)
	{
		for (uint16 i=0; i<Element::nDim(); i++)
			ptr[i] += get_qi_n(n, i);
	}
}
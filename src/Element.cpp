#include "Element.h"

uint16 Element::number_of_dimensions=3;
uint16 Element::number_of_dofs = 0;
uint16 Node::number_of_dofs = 3;
uint16 Element::number_of_nodes = 8;
uint16 Element::number_of_integration_points = 2; 
vector<Face> Element::faces;
vector<Dof_Type> Node::dof_types;
vector<Dof_Type> Element::dof_types;


uint16 Face::num_nodes_per_face = 4;

//-------------Face----------------
std::ostream& operator<< (std::ostream& stream,const Face& obj)
{
	for (uint16 i = 0; i < obj.num_nodes_per_face; i++)
		stream << obj.nodes[i] << " ";
	return stream;
}

//--------------Element---------------
void Element::read_from_stream (istream &str)
{
	assert(nodes);
	for (uint16 i=0; i < number_of_nodes; i++)
		str >> nodes[i];  //TODO: process wrong input
}

void Element::display (uint32 en)
{
	char buffer[256]="";
	ostrstream str(buffer, 256);
	str << "E " << en << ":";
	for (uint16 i=0; i < n_nodes(); i++)
	{
		str << node_num(i);
		if (i+1 != number_of_nodes)
			str << "%t";
	}
	echo(buffer);
}

string Element::toString()
{
	char buffer[256] = "";
	ostrstream str(buffer,256);
	for (uint16 i=0; i < n_nodes(); i++)
	{
		str << node_num(i);
		if (i+1 != number_of_nodes)
			str << " ";
	}
	return string(buffer);
}

Element& Element::operator= (const Element& from)
{
	assert(nodes && from.nodes);
	memcpy(nodes, from.nodes, sizeof(uint32)*number_of_nodes);
	return *this;
}

void Element::change_node_dofs_num (uint16 ndof,...)
{
	Node::number_of_dofs = ndof;
	if (Node::dof_types.size() == 0)
	{
		va_list vlist;
		va_start(vlist, ndof);
		for (uint16 i=0; i < ndof; i++)
			Node::dof_types.push_back(va_arg(vlist,Dof_Type));
	}
}
void Element::change_el_dofs_num (uint16 ndof, ...)
{
	//TODO: должно выполнятся один раз!
	Element::number_of_dofs = ndof;
	if (Element::dof_types.size() == 0)
	{
		va_list vlist;
		va_start(vlist, ndof);
		for (uint16 i=0; i < ndof; i++)
			Element::dof_types.push_back(va_arg(vlist,Dof_Type));
	}
}

Element::elTypes Element::elName2elType (string elName) {
  for (uint16 i = 0; i < Element::LAST; i++) {
    if (elName.compare(Element::elTypeLabels[i]) == 0) {
      return i;
    }
  }
  return Element::NOT_DEFINED;
}
void Element::createElements (elTypes elId, const uint32 n, void** ptr) {
  if (elId == Element::NOT_DEFINED)
    error("createElements: can't find element type with name %s", elName.c_str());
  switch (elId) {
    case Element::PLANE41:
      for (uint32 i = 0; i < n; i++) {
        ptr[i] = (void*) new MIXED_4N_2D_P0();
      }
      break;
    case Element::SOLID81:
      for (uint32 i = 0; i < n; i++) {
        ptr[i] = (void*) new MIXED_8N_3D_P0();
      }
      break;
    default:
      error("createElements: don't have an element with id %d", matId);
  }
}


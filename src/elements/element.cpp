// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "elements/element.h"

namespace nla3d {

FEStorage* Element::storage = NULL;
uint16 Element::number_of_dimensions=3;
uint16 Element::number_of_dofs = 0;
uint16 Node::number_of_dofs = 3;
uint16 Element::number_of_nodes = 8;
uint16 Element::number_of_integration_points = 2; 
std::vector<Face> Element::faces;
std::vector<Dof_Type> Node::dof_types;
std::vector<Dof_Type> Element::dof_types;


uint16 Face::num_nodes_per_face = 4;

//-------------Face----------------
std::ostream& operator<< (std::ostream& stream,const Face& obj)
{
	for (uint16 i = 0; i < obj.num_nodes_per_face; i++)
		stream << obj.nodes[i] << " ";
	return stream;
}

//--------------Element---------------
void Element::read_from_stream (std::istream &str)
{
	assert(nodes);
	for (uint16 i=0; i < number_of_nodes; i++)
		str >> nodes[i];  //TODO: process wrong input
}

void Element::print (std::ostream& out) {
	out << "E " << getElNum() << ":";
	for (uint16 i = 0; i < n_nodes(); i++) {
		out << "\t" << node_num(i);
	}
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
  va_list vlist;
	if (Node::dof_types.size() == 0) {
		va_start(vlist, ndof);
		for (uint16 i=0; i < ndof; i++) {
			Node::dof_types.push_back((Dof_Type)va_arg(vlist, int));
    }
	}
}
void Element::change_el_dofs_num (uint16 ndof, ...)
{
	//TODO: должно выполнятся один раз!
	Element::number_of_dofs = ndof;
  va_list vlist;
	if (Element::dof_types.size() == 0) {
		va_start(vlist, ndof);
		for (uint16 i=0; i < ndof; i++)
			Element::dof_types.push_back((Dof_Type)va_arg(vlist, int));
	}
}

} // namespace nla3d

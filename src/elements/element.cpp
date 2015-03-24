// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d 

#include "elements/element.h"

namespace nla3d {

FEStorage* Element::storage = NULL;
uint16 Element::number_of_dimensions=3;
uint16 Element::number_of_nodes = 8;
uint16 Element::number_of_integration_points = 2; 
std::vector<Face> Element::faces;
std::vector<uint16> Element::dofNumberList;
uint16 Element::numberOfDofs = 0;

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
		out << "\t" << getNodeNumber(i);
	}
}

Element& Element::operator= (const Element& from)
{
	assert(nodes && from.nodes);
	memcpy(nodes, from.nodes, sizeof(uint32)*number_of_nodes);
	return *this;
}

void Element::registerDofType(Dof::dofType type) {
  bool isFound = false;
  if (dofNumberList[type] == Dof::UNDEFINED) {
    dofNumberList[type] = numberOfDofs;
    numberOfDofs++;
  }
}

Dof::dofType Element::getDofType (uint16 dofIndex) {
  assert(dofIndex < getNumberOfDofs());
  for (size_t i = 0; i < Dof::numberOfDofTypes; i++) {
    if (getDofIndex((Dof::dofType)i) == dofIndex) {
      return (Dof::dofType) i;
    }
  }
  error ("What I'm doing here?");
  return Dof::UNDEFINED;
}
} // namespace nla3d

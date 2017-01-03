#include "sys.h"
#include "FEStorage.h"
#include "elements/element.h"
#include "elements/isoparametric.h"

//TODO: not all results are autochecked
static const double th = 1.0e-8; 

using namespace std;
using namespace nla3d;
using namespace nla3d::math;

class LINEdummy : public ElementIsoParamLINE {
    virtual void pre() {
      makeJacob();
      storage->addElementDof(elNum, {Dof::UX});
    };
    virtual void buildK() { };
    virtual void update() { };
};

class QUADdummy : public ElementIsoParamQUAD {
    virtual void pre() {
      makeJacob();
      storage->addElementDof(elNum, {Dof::UX});
    };
    virtual void buildK() { };
    virtual void update() { };
};

class HEXAHEDRONdummy : public ElementIsoParamHEXAHEDRON {
    virtual void pre() {
      makeJacob();
      storage->addElementDof(elNum, {Dof::UX});
    };
    virtual void buildK() { };
    virtual void update() { };
};


int main() {

  cout << "Checking ElementIsoParamLINE..." << endl;
  for (uint16 intOrder = 0; intOrder <= 3; intOrder++) {
    cout << "  Checking for " << intOrder << "number of int. order..." << endl;
    const uint32 numberOfNodes = 8;
    double nodeTable[numberOfNodes][3] = {
// an element like in local coord system
                          {-1.0,  0.0, 0.0},       // 1
                          {+1.0,  0.0, 0.0},       // 2
// an element like in local coord system, but stretched x2
                          {-2.0,  0.0, 0.0},       // 3
                          {+2.0,  0.0, 0.0},       // 4
// arbitrary element nodes
                          {-12.0, 1.0, 3.0},       // 5
                          {4.0, -5.0, -5.0},       // 6
// arbitrary element nodes translated by 
                          {-12.0 + 9.0, +1.0 - 5.0, 3.0 - 3.0},   // 7
                          {+4.0  + 9.0, -5.0 -5.0, -5.0 - 3.0}};  // 8

  // Create an instance of FEStorage.
	FEStorage storage;
  // Create and add nodes into FEStorage
  for (uint32 i = 1; i <= numberOfNodes; i++) {
    Node* no = new Node;
    no->pos[0] = nodeTable[i-1][0];
    no->pos[1] = nodeTable[i-1][1];
    no->pos[2] = nodeTable[i-1][2];
    storage.addNode(no);
  }

  // Create elements instances, define needed element parameters and add them into FEStorage.
  LINEdummy* el = new LINEdummy;
  el->getNodeNumber(0) = 1;
  el->getNodeNumber(1) = 2;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 1

  el = new LINEdummy;
  el->getNodeNumber(0) = 3;
  el->getNodeNumber(1) = 4;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 2

  // arbitrary node position
  el = new LINEdummy;
  el->getNodeNumber(0) = 5;
  el->getNodeNumber(1) = 6;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 3

  // translated node positions
  el = new LINEdummy;
  el->getNodeNumber(0) = 7;
  el->getNodeNumber(1) = 8;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 4

  // switch node numbers
  el = new LINEdummy;
  el->getNodeNumber(0) = 6;  
  el->getNodeNumber(1) = 5;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 5

  // toggle Element::pre() for all
  storage.initializeSolutionData();

  cout << "Checking ElementIsoParamLINE (-1),(+1) element.." << endl;
  el = dynamic_cast<LINEdummy*>(&storage.getElement(1));
  cout << "el->det : " << el->det << endl;
  CHECK_EQTH(el->det, 1.0, th);
  cout << "Volume = " << el->volume() << endl;
  CHECK_EQTH(el->volume(), 2.0, th);


  cout << "Checking ElementIsoParamLINE (-2),(+2) element.." << endl;
  el = dynamic_cast<LINEdummy*>(&storage.getElement(2));
  cout << "el->det : " << el->det << endl;
  CHECK_EQTH(el->det, 2.0, th);
  cout << "Volume = " << el->volume() << endl;
  CHECK_EQTH(el->volume(), 4.0, th);


  cout << "Checking ElementIsoParamLINE with arbitrary nodes .." << endl;
  LINEdummy* el1 = dynamic_cast<LINEdummy*>(&storage.getElement(3));
  LINEdummy* el2 = dynamic_cast<LINEdummy*>(&storage.getElement(4));
  LINEdummy* el3 = dynamic_cast<LINEdummy*>(&storage.getElement(5));
  assert(el1->getIntegrationOrder() == el2->getIntegrationOrder());
  assert(el2->getIntegrationOrder() == el3->getIntegrationOrder());
  cout << "el1->det : " << el1->det << endl;
  cout << "el2->det : " << el2->det << endl;
  cout << "el3->det : " << el3->det << endl;
  CHECK_EQTH(el1->det, el2->det, th);
  CHECK_EQTH(el2->det, el3->det, th);

  cout << "el1->Volume = " << el1->volume() << endl;
  cout << "el2->Volume = " << el2->volume() << endl;
  cout << "el3->Volume = " << el3->volume() << endl;
  CHECK_EQTH(el1->volume(), el2->volume(), th);
  CHECK_EQTH(el2->volume(), el3->volume(), th);
  }


  cout << "Checking ElementIsoParamQUAD..." << endl;
  for (uint16 intOrder = 0; intOrder <= 3; intOrder++) {
    cout << "  Checking for " << intOrder << "number of int. order..." << endl;
    const uint32 numberOfNodes = 16;
    double nodeTable[numberOfNodes][3] = {
// an element like in local coord system
                          {-1.0, -1.0, 0.0},       // 1
                          {1.0, -1.0, 0.0},        // 2
                          {1.0, 1.0, 0.0},         // 3
                          {-1.0, 1.0, 0.0},        // 4
// an element like in local coord system, but stretched x2
                          {-1.0 * 2, -1.0, 0.0},   // 5
                          {1.0 * 2, -1.0, 0.0},    // 6
                          {1.0 * 2, 1.0, 0.0},     // 7
                          {-1.0 * 2, 1.0, 0.0},    // 8
// arbitrary element nodes
                          {-12.0, 1.0, 0.0},       // 9
                          {4.0, -5.0, 0.0},        // 10
                          {7.0, 1.0, 0.0},         // 11
                          {-2.0, 3.0, 0.0},        // 12
// arbitrary element nodes translated by 
                          {-12.0 + 9.0, 1.0 - 5.0, 0.0},       // 13
                          {4.0 + 9.0, -5.0 - 5.0, 0.0},        // 14
                          {7.0 + 9.0, 1.0 - 5.0, 0.0},         // 15
                          {-2.0 + 9.0, 3.0 - 5.0, 0.0}};       // 16

  // Create an instance of FEStorage.
	FEStorage storage;
  // Create and add nodes into FEStorage
  for (uint32 i = 1; i <= numberOfNodes; i++) {
    Node* no = new Node;
    no->pos[0] = nodeTable[i-1][0];
    no->pos[1] = nodeTable[i-1][1];
    no->pos[2] = nodeTable[i-1][2];
    storage.addNode(no);
  }

  // Create elements instances, define needed element parameters and add them into FEStorage.
  QUADdummy* el = new QUADdummy;
  el->getNodeNumber(0) = 1;
  el->getNodeNumber(1) = 2;
  el->getNodeNumber(2) = 3;
  el->getNodeNumber(3) = 4;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 1

  el = new QUADdummy;
  el->getNodeNumber(0) = 5;
  el->getNodeNumber(1) = 6;
  el->getNodeNumber(2) = 7;
  el->getNodeNumber(3) = 8;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 2

  // arbitrary node position
  el = new QUADdummy;
  el->getNodeNumber(0) = 9;
  el->getNodeNumber(1) = 10;
  el->getNodeNumber(2) = 11;
  el->getNodeNumber(3) = 12;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 3

  // translated node positions
  el = new QUADdummy;
  el->getNodeNumber(0) = 13;
  el->getNodeNumber(1) = 14;
  el->getNodeNumber(2) = 15;
  el->getNodeNumber(3) = 16;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 4

  // switch node numbers
  el = new QUADdummy;
  el->getNodeNumber(0) = 10;  
  el->getNodeNumber(1) = 11;
  el->getNodeNumber(2) = 12;
  el->getNodeNumber(3) = 9;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 5

  // toggle Element::pre() for all
  storage.initializeSolutionData();

  cout << "Checking ElementIsoParamQUAD (-1,-1),(1,-1),(1,1),(-1,1) element.." << endl;
  el = dynamic_cast<QUADdummy*>(&storage.getElement(1));
  for (uint32 np = 0; np < el->nOfIntPoints(); np++) { 
    cout << "el->det[" << np << "] : " << el->det[np] << endl;
    CHECK_EQTH(el->det[np], 1.0, th);
  }
  cout << "Volume = " << el->volume() << endl;
  CHECK_EQTH(el->volume(), 4.0, th);


  cout << "Checking ElementIsoParamQUAD (-2,-1),(2,-1),(2,1),(-2,1) element.." << endl;
  el = dynamic_cast<QUADdummy*>(&storage.getElement(2));
  for (uint32 np = 0; np < el->nOfIntPoints(); np++) { 
    cout << "el->det[" << np << "] : " << el->det[np] << endl;
    CHECK_EQTH(el->det[np], 2.0, th);
  }
  cout << "Volume = " << el->volume() << endl;
  CHECK_EQTH(el->volume(), 8.0, th);


  cout << "Checking ElementIsoParamQUAD with arbitrary nodes .." << endl;
  QUADdummy* el1 = dynamic_cast<QUADdummy*>(&storage.getElement(3));
  QUADdummy* el2 = dynamic_cast<QUADdummy*>(&storage.getElement(4));
  QUADdummy* el3 = dynamic_cast<QUADdummy*>(&storage.getElement(5));
  assert(el1->getIntegrationOrder() == el2->getIntegrationOrder());
  assert(el2->getIntegrationOrder() == el3->getIntegrationOrder());
  for (uint32 np = 0; np < el1->nOfIntPoints(); np++) { 
    cout << "el1->det[" << np << "] : " << el1->det[np] << endl;
      cout << "el2->det[" << np << "] : " << el2->det[np] << endl;
      cout << "el3->det[" << np << "] : " << el3->det[np] << endl;
      CHECK_EQTH(el1->det[np], el2->det[np], th);

      for (uint16 no = 0; no < 4; no++) 
        for (uint16 dim = 0; dim < 2; dim++) {
          // cout << "el1->NiXj[" << np << "][" << no << "][" << dim << "] " << el1->NiXj[np][no][dim] << endl;
          // cout << "el2->NiXj[" << np << "][" << no << "][" << dim << "] " << el2->NiXj[np][no][dim] << endl;
          // cout << "el3->NiXj[" << np << "][" << no << "][" << dim << "] " << el3->NiXj[np][no][dim] << endl;
          CHECK_EQTH(el1->NiXj[np][no][dim], el2->NiXj[np][no][dim], th);
        }
    }
      cout << "el1->Volume = " << el1->volume() << endl;
      cout << "el2->Volume = " << el2->volume() << endl;
      cout << "el3->Volume = " << el3->volume() << endl;
      CHECK_EQTH(el1->volume(), el2->volume(), th);
      CHECK_EQTH(el2->volume(), el3->volume(), th);
  }
 



  cout << "Checking ElementIsoParamHEXAHEDRON..." << endl;
  for (uint16 intOrder = 0; intOrder <= 3; intOrder++) {
    cout << "  Checking for " << intOrder << "number of int. order..." << endl;
    const uint32 numberOfNodes = 32;
    double nodeTable[numberOfNodes][3] = {
// an element like in local coord system
                          {-1.0, -1.0, -1.0},       // 1
                          {1.0, -1.0, -1.0},        // 2
                          {1.0, 1.0, -1.0},         // 3
                          {-1.0, 1.0, -1.0},        // 4
                          {-1.0, -1.0, 1.0},        // 5
                          {1.0, -1.0, 1.0},         // 6
                          {1.0, 1.0, 1.0},          // 7
                          {-1.0, 1.0, 1.0},         // 8
// an element like in local coord system, but stretched x2
                          {-1.0 * 2, -1.0, -1.0},   // 9
                          {1.0 * 2, -1.0, -1.0},    // 10
                          {1.0 * 2, 1.0, -1.0},     // 11
                          {-1.0 * 2, 1.0, -1.0},    // 12
                          {-1.0 * 2, -1.0, 1.0},    // 13
                          {1.0 * 2, -1.0, 1.0},     // 14
                          {1.0 * 2, 1.0, 1.0},      // 15
                          {-1.0 * 2, 1.0, 1.0},     // 16
// arbitrary element nodes
                          {-12.0, 1.0, -7.0},       // 17
                          {4.0, -5.0, -7.0},        // 18
                          {7.0, 1.0, -7.0},         // 19
                          {-2.0, 3.0, -7.0},        // 20
                          {-12.0, 1.0, 4.0},        // 21
                          {4.0, -5.0, 4.0},         // 22
                          {7.0, 1.0, 4.0},          // 23
                          {-2.0, 3.0, 4.0},         // 24
// arbitrary element nodes translated by 
                          {-12.0 + 9.0, 1.0 - 5.0, -7.0 + 2.0},       // 25
                          {4.0 + 9.0, -5.0 - 5.0, -7.0 + 2.0},        // 26
                          {7.0 + 9.0, 1.0 - 5.0, -7.0 + 2.0},         // 27
                          {-2.0 + 9.0, 3.0 - 5.0, -7.0 + 2.0},       // 28
                          {-12.0 + 9.0, 1.0 - 5.0, 4.0 + 2.0},        // 29
                          {4.0 + 9.0, -5.0 - 5.0, 4.0 + 2.0},         // 30
                          {7.0 + 9.0, 1.0 - 5.0, 4.0 + 2.0},          // 31
                          {-2.0 + 9.0, 3.0 - 5.0, 4.0 + 2.0}};        // 32

  // Create an instance of FEStorage.
	FEStorage storage;
  // Create and add nodes into FEStorage
  for (uint32 i = 1; i <= numberOfNodes; i++) {
    Node* no = new Node;
    no->pos[0] = nodeTable[i-1][0];
    no->pos[1] = nodeTable[i-1][1];
    no->pos[2] = nodeTable[i-1][2];
    storage.addNode(no);
  }

  // Create elements instances, define needed element parameters and add them into FEStorage.
  HEXAHEDRONdummy* el = new HEXAHEDRONdummy;
  el->getNodeNumber(0) = 1;
  el->getNodeNumber(1) = 2;
  el->getNodeNumber(2) = 3;
  el->getNodeNumber(3) = 4;
  el->getNodeNumber(4) = 5;
  el->getNodeNumber(5) = 6;
  el->getNodeNumber(6) = 7;
  el->getNodeNumber(7) = 8;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 1

  el = new HEXAHEDRONdummy;
  el->getNodeNumber(0) = 9;
  el->getNodeNumber(1) = 10;
  el->getNodeNumber(2) = 11;
  el->getNodeNumber(3) = 12;
  el->getNodeNumber(4) = 13;
  el->getNodeNumber(5) = 14;
  el->getNodeNumber(6) = 15;
  el->getNodeNumber(7) = 16;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 2

  // arbitrary node position
  el = new HEXAHEDRONdummy;
  el->getNodeNumber(0) = 17;
  el->getNodeNumber(1) = 18;
  el->getNodeNumber(2) = 19;
  el->getNodeNumber(3) = 20;
  el->getNodeNumber(4) = 21;
  el->getNodeNumber(5) = 22;
  el->getNodeNumber(6) = 23;
  el->getNodeNumber(7) = 24;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 3

  // translated node positions
  el = new HEXAHEDRONdummy;
  el->getNodeNumber(0) = 25;
  el->getNodeNumber(1) = 26;
  el->getNodeNumber(2) = 27;
  el->getNodeNumber(3) = 28;
  el->getNodeNumber(4) = 29;
  el->getNodeNumber(5) = 30;
  el->getNodeNumber(6) = 31;
  el->getNodeNumber(7) = 32;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 4

  // switch node numbers
  el = new HEXAHEDRONdummy;
  el->getNodeNumber(0) = 20;
  el->getNodeNumber(1) = 17;
  el->getNodeNumber(2) = 18;
  el->getNodeNumber(3) = 19;
  el->getNodeNumber(4) = 24;
  el->getNodeNumber(5) = 21;
  el->getNodeNumber(6) = 22;
  el->getNodeNumber(7) = 23;
  el->setIntegrationOrder(intOrder);
  storage.addElement(el);       // 5

  // toggle Element::pre() for all
  storage.initializeSolutionData();

  cout << "Checking ElementIsoParamHEXAHEDRON uni element.." << endl;
  el = dynamic_cast<HEXAHEDRONdummy*>(&storage.getElement(1));
  for (uint32 np = 0; np < el->nOfIntPoints(); np++) { 
    cout << "el->det[" << np << "] : " << el->det[np] << endl;
    CHECK_EQTH(el->det[np], 1.0, th);
  }
  cout << "Volume = " << el->volume() << endl;
  CHECK_EQTH(el->volume(), 8.0, th);


  cout << "Checking ElementIsoParamHEXAHEDRON stretched uni element.." << endl;
  el = dynamic_cast<HEXAHEDRONdummy*>(&storage.getElement(2));
  for (uint32 np = 0; np < el->nOfIntPoints(); np++) { 
    cout << "el->det[" << np << "] : " << el->det[np] << endl;
    CHECK_EQTH(el->det[np], 2.0, th);
  }
  cout << "Volume = " << el->volume() << endl;
  CHECK_EQTH(el->volume(), 16.0, th);


  cout << "Checking ElementIsoParamHEXAHEDRON with arbitrary nodes .." << endl;
  HEXAHEDRONdummy* el1 = dynamic_cast<HEXAHEDRONdummy*>(&storage.getElement(3));
  HEXAHEDRONdummy* el2 = dynamic_cast<HEXAHEDRONdummy*>(&storage.getElement(4));
  HEXAHEDRONdummy* el3 = dynamic_cast<HEXAHEDRONdummy*>(&storage.getElement(5));
  assert(el1->getIntegrationOrder() == el2->getIntegrationOrder());
  assert(el2->getIntegrationOrder() == el3->getIntegrationOrder());
  for (uint32 np = 0; np < el1->nOfIntPoints(); np++) { 
    cout << "el1->det[" << np << "] : " << el1->det[np] << endl;
      cout << "el2->det[" << np << "] : " << el2->det[np] << endl;
      cout << "el3->det[" << np << "] : " << el3->det[np] << endl;
      CHECK_EQTH(el1->det[np], el2->det[np], th);

      for (uint16 no = 0; no < 8; no++) 
        for (uint16 dim = 0; dim < 3; dim++) {
          // cout << "el1->NiXj[" << np << "][" << no << "][" << dim << "] " << el1->NiXj[np][no][dim] << endl;
          // cout << "el2->NiXj[" << np << "][" << no << "][" << dim << "] " << el2->NiXj[np][no][dim] << endl;
          // cout << "el3->NiXj[" << np << "][" << no << "][" << dim << "] " << el3->NiXj[np][no][dim] << endl;
          CHECK_EQTH(el1->NiXj[np][no][dim], el2->NiXj[np][no][dim], th);
        }
    }
      cout << "el1->Volume = " << el1->volume() << endl;
      cout << "el2->Volume = " << el2->volume() << endl;
      cout << "el3->Volume = " << el3->volume() << endl;
      CHECK_EQTH(el1->volume(), el2->volume(), th);
      CHECK_EQTH(el2->volume(), el3->volume(), th);
    }

}


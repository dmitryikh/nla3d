#include "sys.h"
#include "FESolver.h"

using namespace nla3d;

int main () {
  TimeControl tc;

  double dt = 15.0;
  // testcase:
  // Loadstep 1: dt = 15.0, currentTime = 15.0, numberOfEquilibrium = 1
  // Loadstep 2: dt = 15.0, currentTime = 30.0, numberOfEquilibrium = 2
  // Loadstep 3: dt = 15.0, currentTime = 45.0, numberOfEquilibrium = 3
  // Loadstep 4: dt = 15.0, currentTime = 60.0, numberOfEquilibrium = 4
  // Loadstep 5: dt = 15.0, currentTime = 75.0, numberOfEquilibrium = 5
  // Loadstep 6: dt = 12.5, currentTime = 87.5, numberOfEquilibrium = 6
  // Loadstep 7: dt = 12.5, currentTime = 100 , numberOfEquilibrium = 7
  
  tc.setStartTime(0.0); 
  tc.setEndTime(100.0);
  CHECK(tc.getEndTime() == 100.0);
  CHECK(tc.getStartTime() == 0.0);

  // Loadstep 1
  tc.nextStep(dt);
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 1);
  CHECK(tc.getNumberOfConvergedSteps() == 0);
  CHECK(tc.getCurrentEquilibriumStep() == 1);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 1);
  CHECK(tc.getCurrentTime() == 15.0);
  CHECK(tc.getCurrentNormalizedTime() == 0.15);
  CHECK(tc.getCurrentTimeDelta() == 15.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.15);

  // Loadstep 2
  tc.nextStep(dt);
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 2);
  CHECK(tc.getNumberOfConvergedSteps() == 1);
  CHECK(tc.getCurrentEquilibriumStep() == 1);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 2);
  CHECK(tc.getCurrentTime() == 30.0);
  CHECK(tc.getCurrentNormalizedTime() == 0.30);
  CHECK(tc.getCurrentTimeDelta() == 15.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.15);
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 2);
  CHECK(tc.getNumberOfConvergedSteps() == 1);
  CHECK(tc.getCurrentEquilibriumStep() == 2);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 3);
  CHECK(tc.getCurrentTime() == 30.0);
  CHECK(tc.getCurrentNormalizedTime() == 0.30);
  CHECK(tc.getCurrentTimeDelta() == 0.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.0);

  // Loadstep 3
  tc.nextStep(dt);
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 3);
  CHECK(tc.getNumberOfConvergedSteps() == 2);
  CHECK(tc.getCurrentEquilibriumStep() == 1);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 4);
  CHECK(tc.getCurrentTime() == 45.0);
  CHECK(tc.getCurrentNormalizedTime() == 0.45);
  CHECK(tc.getCurrentTimeDelta() == 15.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.15);
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 3);
  CHECK(tc.getNumberOfConvergedSteps() == 2);
  CHECK(tc.getCurrentEquilibriumStep() == 3);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 6);
  CHECK(tc.getCurrentTime() == 45.0);
  CHECK(tc.getCurrentNormalizedTime() == 0.45);
  CHECK(tc.getCurrentTimeDelta() == 0.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.0);

  // Loadstep 4
  tc.nextStep(dt);
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 4);
  CHECK(tc.getNumberOfConvergedSteps() == 3);
  CHECK(tc.getCurrentEquilibriumStep() == 1);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 7);
  CHECK(tc.getCurrentTime() == 60.0);
  CHECK(tc.getCurrentNormalizedTime() == 0.60);
  CHECK(tc.getCurrentTimeDelta() == 15.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.15);
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 4);
  CHECK(tc.getNumberOfConvergedSteps() == 3);
  CHECK(tc.getCurrentEquilibriumStep() == 4);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 10);
  CHECK(tc.getCurrentTime() == 60.0);
  CHECK(tc.getCurrentNormalizedTime() == 0.60);
  CHECK(tc.getCurrentTimeDelta() == 0.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.0);

  // Loadstep 5
  tc.nextStep(dt);
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 5);
  CHECK(tc.getNumberOfConvergedSteps() == 4);
  CHECK(tc.getCurrentEquilibriumStep() == 1);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 11);
  CHECK(tc.getCurrentTime() == 75.0);
  CHECK(tc.getCurrentNormalizedTime() == 0.75);
  CHECK(tc.getCurrentTimeDelta() == 15.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.15);
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 5);
  CHECK(tc.getNumberOfConvergedSteps() == 4);
  CHECK(tc.getCurrentEquilibriumStep() == 5);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 15);
  CHECK(tc.getCurrentTime() == 75.0);
  CHECK(tc.getCurrentNormalizedTime() == 0.75);
  CHECK(tc.getCurrentTimeDelta() == 0.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.0);

  // Loadstep 6
  tc.nextStep(dt);
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 6);
  CHECK(tc.getNumberOfConvergedSteps() == 5);
  CHECK(tc.getCurrentEquilibriumStep() == 1);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 16);
  CHECK(tc.getCurrentTime() == 87.5);
  CHECK(tc.getCurrentNormalizedTime() == 0.875);
  CHECK(tc.getCurrentTimeDelta() == 12.5);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.125);
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 6);
  CHECK(tc.getNumberOfConvergedSteps() == 5);
  CHECK(tc.getCurrentEquilibriumStep() == 6);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 21);
  CHECK(tc.getCurrentTime() == 87.5);
  CHECK(tc.getCurrentNormalizedTime() == 0.875);
  CHECK(tc.getCurrentTimeDelta() == 0.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.0);

  // Loadstep 7
  tc.nextStep(dt);
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 7);
  CHECK(tc.getNumberOfConvergedSteps() == 6);
  CHECK(tc.getCurrentEquilibriumStep() == 1);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 22);
  CHECK(tc.getCurrentTime() == 100);
  CHECK(tc.getCurrentNormalizedTime() == 1.0);
  CHECK(tc.getCurrentTimeDelta() == 12.5);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.125);
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  tc.nextEquilibriumStep();
  CHECK(tc.getCurrentStep() == 7);
  CHECK(tc.getNumberOfConvergedSteps() == 6);
  CHECK(tc.getCurrentEquilibriumStep() == 7);
  CHECK(tc.getTotalNumberOfEquilibriumSteps() == 28);
  CHECK(tc.getCurrentTime() == 100);
  CHECK(tc.getCurrentNormalizedTime() == 1.0);
  CHECK(tc.getCurrentTimeDelta() == 0.0);
  CHECK(tc.getCurrentNormalizedTimeDelta() == 0.0);

  // next step to no where..
  CHECK(tc.nextStep(dt) == false);
}

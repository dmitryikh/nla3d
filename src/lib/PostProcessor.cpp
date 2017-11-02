// This file is a part of nla3d project. For information about authors and
// licensing go to project's repository on github:
// https://github.com/dmitryikh/nla3d

#include "PostProcessor.h"

namespace nla3d {

PostProcessor::PostProcessor(FEStorage *st) {
  assert(st);
  storage = st;
  active = false;
  failed = false;
}

std::string PostProcessor::getStatus () {
  char buf[300];
  sprintf_s(buf, 200, "PostProcessor No %d: %s\n\tActive - %s, Failed - %s",nPost_proc, name.c_str(),active?"true":"false",active?"true":"false");
  return std::string(buf);
}

} // namespace nla3d

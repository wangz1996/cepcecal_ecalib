#include "TreeManager.hh"

int main(int argc, char** argv) {
    auto tm = new TreeManager();

    tm->readTree("/mnt2/USTC/jxwang/CEPC_ScECAL/CEPCforZhen/Calib_IncludeADC/40GeV_hl_electron.root","Calib_Hit");
    tm->eventLoop();

    return 1;
}
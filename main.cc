#include "TreeManager.hh"

int main(int argc, char** argv) {
    auto tm = new TreeManager();
    for(size_t i=0; i<argc; i++) {
        if(strcmp(argv[i],"-s")==0) {
            tm->setDoScale(argv[i+1]);
            std::cout<<"Using scale factor "<<argv[i+1]<<std::endl;
        }
    }
    tm->readTree("/mnt2/USTC/jxwang/CEPC_ScECAL/CEPCforZhen/Calib_IncludeADC/40GeV_hl_electron.root","Calib_Hit");
    tm->eventLoop();

    return 1;
}
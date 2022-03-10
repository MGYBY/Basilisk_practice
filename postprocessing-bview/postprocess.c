#include "run.h"
#include "view.h"
int main(){
  run();
}
event init(t=0){
  double timebgn=0.0;
  double timestp=0.01;
  double timeend=1;
  char namefile[100];
  int sdfs;
  double area,timeload;
  char nameBview[100],textBview[100];
  scalar varLeft[],varRight[],omega[];
  for(timeload=timebgn;timeload<=timeend;timeload+=timestp){
    sprintf(namefile,"dumpfile%.2f",timeload);
    printf("load the file %s!\n",namefile);
    restore(file=namefile);
    sdfs=0;
    foreach(){
      sdfs+=1;
    }
    printf("number of cells are %d.\n",sdfs);
    sprintf(nameBview,"out-bview-level-vorticity-%.2f.png",timeload);
    view(width=1000,height=500);
    clear();
    draw_vof("f",lw=5);
    box();
    cells();
    save(nameBview);
  }
}
event end(t=0.0){
}

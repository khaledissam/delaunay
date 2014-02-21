#include "delaunay.hpp"
#include <iostream>
#include <cmath>

int main(int argc,char*argv[])
{
  
  delaunay d;
  d.add(0,0);
  d.add(1,1);
  d.add(5,0);
  
  if(d.validate())cout<<"The delaunay triangulation is valid\n";
  else cout<<"The delaunay triangulation is not valid\n";
  return 0;
}

//
//  main.cpp
//  myproject
//
//  Created by Jaekwang Kim on 2/13/20.
//  Copyright © 2020 Jaekwang Kim. All rights reserved.
//

#include <iostream>
#include "KWCJumpFunction.h"

using namespace std;

/* this code fits the jump function to match the alternate KWC energy
to given grain boundary energy data, FCC [110] symmetric tilt gradient energy
*/

 int main(int argc, const char * argv[]) {
 
     std::cout << "Program_started"<<std::endl;

     int nData=361;
     /* KWCDesignJ optimizer(number of data,"inputFileName.txt","outputFileName.txt") */;
     KWCDesignJ optimizer(nData,"STGB_cu_110.txt","J_cu_110.txt");
     optimizer.run(); 
   
     return 0;
 }



//
// Created by smdupor on 9/21/21.
//
#include <iostream>
#include <fstream>

int main() {
   std::ofstream fp("dataset.txt");
   std::string temp;

   for(uint64_t i = 0; i < 0b1111111111111111111111; i++){
      temp = std::to_string(rand() % 0xfffff) +  "\t"  + std::to_string(rand() % 0xfffff) + "\n";
      fp << temp;
   }
   fp.close();
}

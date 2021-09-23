//
// Created by smdupor on 9/21/21.
//
#include <iostream>
#include <fstream>

int main() {
   std::ofstream fp("dataset.txt");
   std::string temp;

   for(uint64_t i = 0; i < 0b11111111111111111111111111; i++){
      temp = std::to_string(rand() % 0xfff) +  "\t"  + std::to_string(rand() % 0xffffff) + "\n";
      fp << temp;
   }
   fp.close();
}

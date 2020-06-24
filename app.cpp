#include <iostream>
#include "profile.hpp"

int main() {
  Profile sam("Sam Drakkila", 30, "New York", "USA", "he/him");

  std::cout << sam.view_profile() << "\n" ;
  std::cout << "\n";
  sam.add_hobby("Masturbating");
  std::cout << sam.view_profile() << "\n";

}

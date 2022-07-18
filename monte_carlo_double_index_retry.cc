#pragma once

#include <unordered_map>

#include "sphere_helpers.cc"

struct Bin {
  int x, y;
  size_t binLength;
  bool inBounds;
  
  void checkInBounds() {
    bool outOfBounds = x < 0;
    outOfBounds |= y < 0;
    outOfBounds |= x >= binLength;
    outOfBounds |= y >= binLength;
    
    inBounds = !outOfBounds;
  }
  
  Bin(int x_, int y_, size_t binLength_):
    x{x_}
  , y{y_}
  , binLength{binLength_}  {
    checkInBounds();
  }
  
  Bin operator+(const Bin& other) {
    return Bin(x + other.x, y + other.y, binLength);
  }
  
  bool operator==(const Bin& other) {
    return x == other.x && y == other.y;
  }
};

struct BinHash {
  size_t operator() (const Bin& bin) const {
    return std::hash<int>()(bin.x) ^ std::hash<int>()(bin.y);
  }
};

using BinPoints = std::unordered_map<Bin, std::vector<size_t>>;

std::unordered_map<std::string, std::string> test1;
std::unordered_map<Bin, std::string, BinHash> test2;
std::unordered_map<Bin, std::vector<size_t>, BinHash> test3;

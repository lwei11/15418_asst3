/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

#include "wireroute.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <vector>

#include <unistd.h>
#include <omp.h>
#include <time.h>

void print_stats(const std::vector<std::vector<int>>& occupancy) {
  int max_occupancy = 0;
  long long total_cost = 0;

  for (const auto& row : occupancy) {
    for (const int count : row) {
      max_occupancy = std::max(max_occupancy, count);
      total_cost += count * count;
    }
  }

  std::cout << "Max occupancy: " << max_occupancy << '\n';
  std::cout << "Total cost: " << total_cost << '\n';
}

void write_output(const std::vector<Wire>& wires, const int num_wires, const std::vector<std::vector<int>>& occupancy, const int dim_x, const int dim_y, const int num_threads, std::string input_filename) {
  if (std::size(input_filename) >= 4 && input_filename.substr(std::size(input_filename) - 4) == ".txt") {
    input_filename.resize(std::size(input_filename) - 4);
  }

  const std::string occupancy_filename = input_filename + "_occupancy_" + std::to_string(num_threads) + ".txt";
  const std::string wires_filename = input_filename + "_wires_" + std::to_string(num_threads) + ".txt";

  std::ofstream out_occupancy(occupancy_filename, std::fstream::out);
  if (!out_occupancy) {
    std::cerr << "Unable to open file: " << occupancy_filename << '\n';
    exit(EXIT_FAILURE);
  }

  out_occupancy << dim_x << ' ' << dim_y << '\n';
  for (const auto& row : occupancy) {
    for (const int count : row) {
      out_occupancy << count << ' ';
    }
    out_occupancy << '\n';
  }

  out_occupancy.close();

  std::ofstream out_wires(wires_filename, std::fstream:: out);
  if (!out_wires) {
    std::cerr << "Unable to open file: " << wires_filename << '\n';
    exit(EXIT_FAILURE);
  }

  out_wires << dim_x << ' ' << dim_y << '\n' << num_wires << '\n';

  for (const auto& [start_x, start_y, end_x, end_y, bend1_x, bend1_y] : wires) {
    out_wires << start_x << ' ' << start_y << ' ' << bend1_x << ' ' << bend1_y << ' ';

    if (start_y == bend1_y) {
    // first bend was horizontal

      if (end_x != bend1_x) {
        // two bends

        out_wires << bend1_x << ' ' << end_y << ' ';
      }
    } else if (start_x == bend1_x) {
      // first bend was vertical

      if(end_y != bend1_y) {
        // two bends

        out_wires << end_x << ' ' << bend1_y << ' ';
      }
    }
    out_wires << end_x << ' ' << end_y << '\n';
  }

  out_wires.close();
}

void compute_occupancy(std::vector<Wire>& wires, std::vector<std::vector<int>>& occupancy) {
  int counter = 0;
  for(auto& wire: wires) {
    // printf("Wire # %d\n", counter);
    // printf("Bend is at (%d, %d)\n", wire.bend1_x, wire.bend1_y);
    int x0 = wire.start_x;
    int x1 = wire.end_x;
    int y0 = wire.start_y;
    int y1 = wire.end_y;
    int sX = (x0 > x1) ? x1 : x0;
    int eX = (x0 > x1) ? x0 : x1;
    bool xflipped = (x0 > x1);
    int sY = (y0 > y1) ? y1 : y0;
    int eY = (y0 > y1) ? y0 : y1;
    bool yflipped = (y0 > y1);
    if( wire.bend1_x == wire.start_x) {  // bend vertically (first go along y axis)
      for(int i=sX; i<=eX; i++) { 
        // update all x-axis (1 horizontal line segment)
        // printf("1 : Updated (%d, %d)\n", i, wire.bend1_y);
        occupancy[wire.bend1_y][i] += 1;
      }
      for(int j=sY; j<=eY; j++) {
        // update all y-axis (2 vertical line segments)
        if (j<wire.bend1_y && (xflipped != yflipped)) {
          // printf("2 : Updated (%d, %d)\n", eX, j);
          occupancy[j][eX] += 1;  // line segment before the bend
        }
        else if (j<wire.bend1_y && (xflipped == yflipped)) {
          // printf("2 : Updated (%d, %d)\n", sX, j);
          occupancy[j][sX] += 1;  // line segment before the bend
        }
        if (j>wire.bend1_y && (xflipped != yflipped)) {
          // printf("3 : Updated (%d, %d)\n", sX, j);
          occupancy[j][sX] += 1;  // line segment after the bend
        } 
        else if (j>wire.bend1_y && (xflipped == yflipped)) {
          // printf("3 : Updated (%d, %d)\n", eX, j);
          occupancy[j][eX] += 1;  // line segment after the bend
        }
      }
    } // bend vertically (first go along y axis)
    else {
      for(int i=sX; i<=eX; i++) { // bend horizontally (first go along x axis)
        // update all x-axis (2 horizontal line segments)
        if (i<wire.bend1_x && (xflipped != yflipped)) {
          // printf("4 : Updated (%d, %d)\n", i, sY);
          occupancy[eY][i] += 1; // line segment before the bend
        } 
        else if (i<wire.bend1_x && (xflipped == yflipped)) {
          // printf("4 : Updated (%d, %d)\n", i, sY);
          occupancy[sY][i] += 1; // line segment before the bend
        }
        if (i>wire.bend1_x && (xflipped != yflipped)) {
          // printf("5 : Updated (%d, %d)\n", i, eY);
          occupancy[sY][i] += 1;  // line segment after the bend
        } 
        else if (i>wire.bend1_x && (xflipped == yflipped)) {
          // printf("5 : Updated (%d, %d)\n", i, eY);
          occupancy[eY][i] += 1;  // line segment after the bend
        }
      }
      for(int j=sY; j<=eY; j++) {
        // update all y-axis (1 vertical line segment)
        // printf("6 : Updated (%d, %d)\n", wire.bend1_x, j);
        occupancy[j][wire.bend1_x] += 1;
      }
    } // bend horizontally (first go along x axis)
    counter++;
    // printf("---------------\n");
  }
}

void init_random(Wire& currWire) {
  int x0 = currWire.start_x;
  int x1 = currWire.end_x;
  int y0 = currWire.start_y;
  int y1 = currWire.end_y;
  srand(std::chrono::steady_clock::now().time_since_epoch().count());
  if (x0 == x1 || y0 == y1) 
    return;
  int xdist = abs(x0 - x1);
  int ydist = abs(y0 - y1);
  if (rand() % 2 == 1) {
    int bend_loc = rand() % xdist + 1 + ((x0 > x1) ? x1 : x0);
    currWire.bend1_x = (x0 > x1)? bend_loc -1 : bend_loc;
    currWire.bend1_y = y0;
    // printf("Entered X: X0 is %d, X1 is %d, Y0 is %d, Y1 is %d, bend1 is (%d, %d)\n", x0, x1, y0, y1, currWire.bend1_x, currWire.bend1_y);
  }
  else {
    int bend_loc = rand() % ydist + 1 + ((y0 > y1) ? y1 : y0);
    currWire.bend1_y = (y0 > y1)? bend_loc -1 : bend_loc;
    currWire.bend1_x = x0;
    // printf("ENTERED Y: X0 is %d, X1 is %d, Y0 is %d, Y1 is %d, bend1 is (%d, %d)\n", x0, x1, y0, y1, currWire.bend1_x, currWire.bend1_y);
  }
}

void decrement_occupancy(Wire& wire, std::vector<std::vector<int>>& occupancy) {
  int x0 = wire.start_x;
  int x1 = wire.end_x;
  int y0 = wire.start_y;
  int y1 = wire.end_y;
  int sX = (x0 > x1) ? x1 : x0;
  int eX = (x0 > x1) ? x0 : x1;
  bool xflipped = (x0 > x1);
  int sY = (y0 > y1) ? y1 : y0;
  int eY = (y0 > y1) ? y0 : y1;
  bool yflipped = (y0 > y1);
  if( wire.bend1_x == wire.start_x) {
    for(int i=sX; i<=eX; i++) { 
      occupancy[wire.bend1_y][i] -= 1;
    }
    for(int j=sY; j<=eY; j++) {
      if (j<wire.bend1_y && (xflipped != yflipped)) {
        occupancy[j][eX] -= 1;
      }
      else if (j<wire.bend1_y && (xflipped == yflipped)) {
        occupancy[j][sX] -= 1;
      }
      if (j>wire.bend1_y && (xflipped != yflipped)) {
        occupancy[j][sX] -= 1;
      } 
      else if (j>wire.bend1_y && (xflipped == yflipped)) {
        occupancy[j][eX] -= 1;
      }
    }
  }
  else {
    for(int i=sX; i<=eX; i++) {
      if (i<wire.bend1_x && (xflipped != yflipped)) {
        occupancy[eY][i] -= 1;
      } 
      else if (i<wire.bend1_x && (xflipped == yflipped)) {
        occupancy[sY][i] -= 1;
      }
      if (i>wire.bend1_x && (xflipped != yflipped)) {
        occupancy[sY][i] -= 1;
      } 
      else if (i>wire.bend1_x && (xflipped == yflipped)) {
        occupancy[eY][i] -= 1;
      }
    }
    for(int j=sY; j<=eY; j++) {
      occupancy[j][wire.bend1_x] -= 1;
    }
  }
}

void increment_occupancy(Wire& wire, std::vector<std::vector<int>>& occupancy) {
  int x0 = wire.start_x;
  int x1 = wire.end_x;
  int y0 = wire.start_y;
  int y1 = wire.end_y;
  int sX = (x0 > x1) ? x1 : x0;
  int eX = (x0 > x1) ? x0 : x1;
  bool xflipped = (x0 > x1);
  int sY = (y0 > y1) ? y1 : y0;
  int eY = (y0 > y1) ? y0 : y1;
  bool yflipped = (y0 > y1);
  if( wire.bend1_x == wire.start_x) {
    for(int i=sX; i<=eX; i++) { 
      occupancy[wire.bend1_y][i] += 1;
    }
    for(int j=sY; j<=eY; j++) {
      if (j<wire.bend1_y && (xflipped != yflipped)) {
        occupancy[j][eX] += 1;
      }
      else if (j<wire.bend1_y && (xflipped == yflipped)) {
        occupancy[j][sX] += 1;
      }
      if (j>wire.bend1_y && (xflipped != yflipped)) {
        occupancy[j][sX] += 1;
      } 
      else if (j>wire.bend1_y && (xflipped == yflipped)) {
        occupancy[j][eX] += 1;
      }
    }
  }
  else {
    for(int i=sX; i<=eX; i++) {
      if (i<wire.bend1_x && (xflipped != yflipped)) {
        occupancy[eY][i] += 1;
      } 
      else if (i<wire.bend1_x && (xflipped == yflipped)) {
        occupancy[sY][i] += 1;
      }
      if (i>wire.bend1_x && (xflipped != yflipped)) {
        occupancy[sY][i] += 1;
      } 
      else if (i>wire.bend1_x && (xflipped == yflipped)) {
        occupancy[eY][i] += 1;
      }
    }
    for(int j=sY; j<=eY; j++) {
      occupancy[j][wire.bend1_x] += 1;
    }
  }
}

int compute_cost(int x0, int x1, int y0, int y1, int bx, int by, std::vector<std::vector<int>>& occupancy) {
  int cost = 0;
  int sX = (x0 > x1) ? x1 : x0;
  int eX = (x0 > x1) ? x0 : x1;
  bool xflipped = (x0 > x1);
  int sY = (y0 > y1) ? y1 : y0;
  int eY = (y0 > y1) ? y0 : y1;
  bool yflipped = (y0 > y1);
  if(bx == x0) {
    for(int i=sX; i<=eX; i++) { 
      cost += occupancy[by][i];
    }
    for(int j=sY; j<=eY; j++) {
      if (j<by && (xflipped != yflipped)) {
        cost += occupancy[j][eX];
      }
      else if (j<by && (xflipped == yflipped)) {
        cost += occupancy[j][sX];
      }
      if (j>by && (xflipped != yflipped)) {
        cost += occupancy[j][sX];
      } 
      else if (j>by && (xflipped == yflipped)) {
        cost += occupancy[j][eX];
      }
    }
  }
  else {
    for(int i=sX; i<=eX; i++) {
      if (i<bx && (xflipped != yflipped)) {
        cost += occupancy[eY][i];
      } 
      else if (i<bx && (xflipped == yflipped)) {
        cost += occupancy[sY][i];
      }
      if (i>bx && (xflipped != yflipped)) {
        cost += occupancy[sY][i];
      } 
      else if (i>bx && (xflipped == yflipped)) {
        cost += occupancy[eY][i];
      }
    }
    for(int j=sY; j<=eY; j++) {
      cost += occupancy[j][bx];
    }
  }
  return cost;
}

long long compute_potential_cost(int x0, int x1, int y0, int y1, int bX, int bY, std::vector<std::vector<int>>& occupancy) {
  long long cost = 0;
  int sX = (x0 > x1) ? x1 : x0;
  int eX = (x0 > x1) ? x0 : x1;
  bool xflipped = (x0 > x1);
  int sY = (y0 > y1) ? y1 : y0;
  int eY = (y0 > y1) ? y0 : y1;
  bool yflipped = (y0 > y1);
  if(bX == x0) {
    for(int i=sX; i<=eX; i++) { 
      int temp = occupancy[bY][i];
      cost += temp;
    }
    for(int j=sY; j<=eY; j++) {
      if (j<bY && (xflipped != yflipped)) {
        int temp = occupancy[j][eX];
        cost += temp;
      }
      else if (j<bY && (xflipped == yflipped)) {
        int temp = occupancy[j][sX];
        cost += temp;
      }
      if (j>bY && (xflipped != yflipped)) {
        int temp = occupancy[j][sX];
        cost += temp;
      } 
      else if (j>bY && (xflipped == yflipped)) {
        int temp = occupancy[j][eX];
        cost += temp;
      }
    }
  }
  else {
    for(int i=sX; i<=eX; i++) {
      if (i<bX && (xflipped != yflipped)) {
        int temp = occupancy[eY][i];
        cost += temp;
      } 
      else if (i<bX && (xflipped == yflipped)) {
        int temp = occupancy[sY][i];
        cost += temp;
      }
      if (i>bX && (xflipped != yflipped)) {
        int temp = occupancy[sY][i];
        cost += temp;
      } 
      else if (i>bX && (xflipped == yflipped)) {
        int temp = occupancy[eY][i];
        cost += temp;
      }
    }
    for(int j=sY; j<=eY; j++) {
      int temp = occupancy[j][bX];
      cost += temp;
    }
  }
  return cost;
}

void compute_shortest_path(Wire& wire, std::vector<std::vector<int>>& occupancy, long long currCost, const double SA_prob) {
  int x0 = wire.start_x;
  int x1 = wire.end_x;
  int y0 = wire.start_y;
  int y1 = wire.end_y;
  int rX = wire.bend1_x;
  int rY = wire.bend1_y;
  int sX = (x0 > x1) ? x1 : x0;
  int eX = (x0 > x1) ? x0 : x1;
  bool xflipped = (x0 > x1);
  int sY = (y0 > y1) ? y1 : y0;
  int eY = (y0 > y1) ? y0 : y1;
  bool yflipped = (y0 > y1);
  long long finCost = currCost;

  //Decrement original matrix
  if( rX == x0) {
    for(int i=sX; i<=eX; i++) { 
      occupancy[rY][i] -= 1;
    }
    for(int j=sY; j<=eY; j++) {
      if (j<rY && (xflipped != yflipped)) {
        occupancy[j][eX] -= 1;
      }
      else if (j<rY && (xflipped == yflipped)) {
        occupancy[j][sX] -= 1;
      }
      if (j>rY && (xflipped != yflipped)) {
        occupancy[j][sX] -= 1;
      } 
      else if (j>rY && (xflipped == yflipped)) {
        occupancy[j][eX] -= 1;
      }
    }
  }
  else {
    for(int i=sX; i<=eX; i++) {
      if (i<rX && (xflipped != yflipped)) {
        occupancy[eY][i] -= 1;
      } 
      else if (i<rX && (xflipped == yflipped)) {
        occupancy[sY][i] -= 1;
      }
      if (i>rX && (xflipped != yflipped)) {
        occupancy[sY][i] -= 1;
      } 
      else if (i>rX && (xflipped == yflipped)) {
        occupancy[eY][i] -= 1;
      }
    }
    for(int j=sY; j<=eY; j++) {
      occupancy[j][rX] -= 1;
    }
  }

  if (rand() % 100 < (SA_prob * 100)) {
    init_random(wire);
    rX = wire.bend1_x;
    rY = wire.bend1_y;
  }
  else {
    long long tempCost;
    
    //Shortest horizontal path
    #pragma omp parallel for shared(finCost, rX, rY)
      for (int i=1;i<=eX-sX;i++) {
        int bendx = xflipped ? (eX - i) : (sX + i);
        tempCost = compute_potential_cost(x0, x1, y0, y1, bendx, y0, occupancy);
        if (tempCost < finCost) {
          #pragma omp critical
          {
            finCost = tempCost;
            rX = bendx;
            rY = y0;
          }
        }
      }

    #pragma omp parallel for shared(finCost, rX, rY)
      //Shortest vertical path
      for (int j=1;j<=eY-sY;j++) {
        int bendy = yflipped ? (eY - j) : (sY + j);
        tempCost = compute_potential_cost(x0, x1, y0, y1, x0, bendy, occupancy);
        if (tempCost < finCost) {
          #pragma omp critical
          {
            finCost = tempCost;
            rX = x0;
            rY = bendy;
          }
        }
      }
    
    //Update Wire and Matrix
  } 

  wire.bend1_x = rX;
  wire.bend1_y = rY;
  if( rX == x0) {  // bend vertically (first go along y axis)
    for(int i=sX; i<=eX; i++) { 
      occupancy[rY][i] += 1;
    }
    for(int j=sY; j<=eY; j++) {
      // update all y-axis (2 vertical line segments)
      if (j<rY && (xflipped != yflipped)) {
        occupancy[j][eX] += 1;  // line segment before the bend
      }
      else if (j<rY && (xflipped == yflipped)) {
        occupancy[j][sX] += 1;  // line segment before the bend
      }
      if (j>rY && (xflipped != yflipped)) {
        occupancy[j][sX] += 1;  // line segment after the bend
      } 
      else if (j>rY && (xflipped == yflipped)) {
        occupancy[j][eX] += 1;
      }
    }
  }
  else {
    for(int i=sX; i<=eX; i++) {
      if (i<rX && (xflipped != yflipped)) {
        occupancy[eY][i] += 1;
      } 
      else if (i<rX && (xflipped == yflipped)) {
        occupancy[sY][i] += 1;
      }
      if (i>rX && (xflipped != yflipped)) {
        occupancy[sY][i] += 1;
      } 
      else if (i>rX && (xflipped == yflipped)) {
        occupancy[eY][i] += 1;
      }
    }
    for(int j=sY; j<=eY; j++) {
      occupancy[j][rX] += 1;
    }
  }
}

void compute_shortest_path_across(Wire& wire, std::vector<std::vector<int>>& occupancy, long long currCost, const double SA_prob) {
  int x0 = wire.start_x;
  int x1 = wire.end_x;
  int y0 = wire.start_y;
  int y1 = wire.end_y;
  int rX = wire.bend1_x;
  int rY = wire.bend1_y;
  int sX = (x0 > x1) ? x1 : x0;
  int eX = (x0 > x1) ? x0 : x1;
  bool xflipped = (x0 > x1);
  int sY = (y0 > y1) ? y1 : y0;
  int eY = (y0 > y1) ? y0 : y1;
  bool yflipped = (y0 > y1);
  long long finCost = currCost;

  #pragma omp critical 
  {
    decrement_occupancy(wire, occupancy);
  }

  if (rand() % 100 < (SA_prob * 100)) {
    init_random(wire);
    rX = wire.bend1_x;
    rY = wire.bend1_y;
  }
  else {
    long long tempCost;
    
    //Shortest horizontal path
    for (int i=1;i<=eX-sX;i++) {
      int bendx = xflipped ? (eX - i) : (sX + i);
      tempCost = compute_potential_cost(x0, x1, y0, y1, bendx, y0, occupancy);
      if (tempCost < finCost) {
        finCost = tempCost;
        rX = bendx;
        rY = y0;
      }
    }

    //Shortest vertical path
    for (int j=1;j<=eY-sY;j++) {
      int bendy = yflipped ? (eY - j) : (sY + j);
      tempCost = compute_potential_cost(x0, x1, y0, y1, x0, bendy, occupancy);
      if (tempCost < finCost) {
        finCost = tempCost;
        rX = x0;
        rY = bendy;
      }
    }
    
    //Update Wire and Matrix
  } 

  wire.bend1_x = rX;
  wire.bend1_y = rY;
  #pragma omp critical 
  {
    increment_occupancy(wire, occupancy);
  }
}

void within_wire(const int num_threads, const double SA_prob, const int SA_iters, 
std::vector<Wire>& wires, std::vector<std::vector<int>>& occupancy)
{
  // printf("currCost is %lld\n", currCost);
  for(int curr_iters = 0; curr_iters < SA_iters; curr_iters ++) {
    for(auto& wire: wires) {
      if(wire.start_x == wire.end_x || wire.start_y == wire.end_y)
        continue; // two end points are in the same row or column
      long long currCost = compute_cost(wire.start_x, wire.end_x, wire.start_y, 
          	    wire.end_y, wire.bend1_x, wire.bend1_y, occupancy);
      compute_shortest_path(wire, occupancy, currCost, SA_prob);
    } // for each wire
  } // SA iterations 
}

void across_wire(const int num_threads, const double SA_prob, const int SA_iters, 
const int batch_size, const int num_wires, std::vector<Wire>& wires, const int dim_x, const int dim_y, std::vector<std::vector<int>>& occupancy) 
{
  int curr_batch_size = batch_size;
  // printf("currCost is %lld\n", currCost);
  for(int curr_iters = 0; curr_iters < SA_iters; curr_iters ++) {
    for(int batch = 0; batch < num_wires / batch_size; batch++) {
      int startWire = batch*batch_size;
      //printf("Entered: Start wire is %d\n", startWire);
      int wireIndex = 0;
      curr_batch_size = (batch_size < num_wires - startWire)? batch_size : num_wires -startWire;
      
      #pragma omp parallel for shared(wires, occupancy) schedule(dynamic)
        for (wireIndex = 0; wireIndex < curr_batch_size; wireIndex++) {
          //printf("Entered wire calc for wire %d\n", wireIndex);
          Wire& wire = wires[wireIndex + startWire];
          if(wire.start_x == wire.end_x || wire.start_y == wire.end_y)
            continue; // two end points are in the same row or column
          long long currCost = compute_cost(wire.start_x, wire.end_x, wire.start_y, wire.end_y, wire.bend1_x, wire.bend1_y, occupancy);
          compute_shortest_path_across(wire, occupancy, currCost, SA_prob);
        }
    } // for each wire
  } // SA iterations 
}

int main(int argc, char *argv[]) {
  const auto init_start = std::chrono::steady_clock::now();

  std::string input_filename;
  int num_threads = 0;
  double SA_prob = 0.1;
  int SA_iters = 5;
  char parallel_mode = '\0';
  int batch_size = 1;

  int opt;
  while ((opt = getopt(argc, argv, "f:n:p:i:m:b:")) != -1) {
    switch (opt) {
      case 'f':
        input_filename = optarg;
        break;
      case 'n':
        num_threads = atoi(optarg);
        break;
      case 'p':
        SA_prob = atof(optarg);
        break;
      case 'i':
        SA_iters = atoi(optarg);
        break;
      case 'm':
        parallel_mode = *optarg;
        break;
      case 'b':
        batch_size = atoi(optarg);
        break;
      default:
        std::cerr << "Usage: " << argv[0] << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
        exit(EXIT_FAILURE);
    }
  }

  // Check if required options are provided
  if (empty(input_filename) || num_threads <= 0 || SA_iters <= 0 || (parallel_mode != 'A' && parallel_mode != 'W') || batch_size <= 0) {
    std::cerr << "Usage: " << argv[0] << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
    exit(EXIT_FAILURE);
  }

  std::cout << "Number of threads: " << num_threads << '\n';
  std::cout << "Simulated annealing probability parameter: " << SA_prob << '\n';
  std::cout << "Simulated annealing iterations: " << SA_iters << '\n';
  std::cout << "Input file: " << input_filename << '\n';
  std::cout << "Parallel mode: " << parallel_mode << '\n';
  std::cout << "Batch size: " << batch_size << '\n';

  std::ifstream fin(input_filename);

  if (!fin) {
    std::cerr << "Unable to open file: " << input_filename << ".\n";
    exit(EXIT_FAILURE);
  }

  int dim_x, dim_y;
  int num_wires;

  /* Read the grid dimension and wire information from file */
  fin >> dim_x >> dim_y >> num_wires;

  std::vector<Wire> wires(num_wires);
  std::vector occupancy(dim_y, std::vector<int>(dim_x));

  for (auto& wire : wires) {
    fin >> wire.start_x >> wire.start_y >> wire.end_x >> wire.end_y;
    wire.bend1_x = wire.start_x;
    wire.bend1_y = wire.start_y;
  }

  /* Initialize any additional data structures needed in the algorithm */

  omp_set_num_threads(num_threads);

  const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
  std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';

  const auto compute_start = std::chrono::steady_clock::now();

  /** 
   * Implement the wire routing algorithm here
   * Feel free to structure the algorithm into different functions
   * Don't use global variables.
   * Use OpenMP to parallelize the algorithm. 
   */
  for(auto& wire: wires) {
    init_random(wire);
  }
  compute_occupancy(wires, occupancy);

  if (parallel_mode == 'W')
    within_wire(num_threads, SA_prob, SA_iters, wires, occupancy);
  else
    across_wire(num_threads, SA_prob, SA_iters, batch_size, num_wires, wires, dim_x, dim_y, occupancy);

  const double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - compute_start).count();
  std::cout << "Computation time (sec): " << compute_time << '\n';

  /* Write wires and occupancy matrix to files */
  print_stats(occupancy);
  write_output(wires, num_wires, occupancy, dim_x, dim_y, num_threads, input_filename);
}

validate_wire_t Wire::to_validate_format(void) const 
{
  /* TODO(student): Implement this if you want to use the wr_checker. */
  /* See wireroute.h for details on validate_wire_t. */
  
  validate_wire_t vw;
  vw.cleanup(); 
  
  vw.num_pts = 0;
  vw.p[0].x = start_x;
  vw.p[0].y = start_y;
  vw.num_pts ++;
  vw.p[vw.num_pts].x = bend1_x;
  vw.p[vw.num_pts].y = bend1_y;
  vw.num_pts ++;
  vw.p[vw.num_pts].x = end_x;
  vw.p[vw.num_pts].y = end_y;

  printf("wire: start %d, %d; bend1 %d, %d; end %d, %d\n", 
    vw.p[0].x, vw.p[0].y, vw.p[1].x, vw.p[1].y, vw.p[2].x, vw.p[2].y);

  return vw;
}
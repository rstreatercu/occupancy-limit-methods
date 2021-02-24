#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<random>
#include<chrono>
#include<string>

using namespace std;

// Define loops, vectors
#define iLOOP for (int i = 0; i < xlink_pos_.size(); ++i) 
#define jLOOP for (int j = 0; j < site_pos_.size(); ++j) 
#define trialLOOP for (int trial = 0; trial < n; ++trial) 
#define VDD vector<vector<double> >
#define VDDD vector<vector<vector<double> > >
#define VD vector<double>
#define VBB vector<vector<bool> >
#define VB vector<bool>
#define VII vector<vector<int> >
#define VI vector<int>
#define TOLERANCE 1e-8

class Simulation {
  private:
    // Physical constants
    const double kon_ = 2000.0; // on-rate
    const double dt_; // timestep
    const double k_ = 10.0; // spring constant
    const double E_0_ = 1.0; // binding energy
    const double tmax_ = 0.2; // end time

    // Useful variables
    int collision_count_ = 0; // how many collisions occur
    VI xlink_index_ = {0,1}; // order in which xlinks interact
    VI site_index_ = {0,1,2}; // order in which sites are interacted with
    ofstream ofile_; // file to write to

    // xlink, site positions as list of (x,y)
    VDD xlink_pos_ = VDD(2, VD(2,0.0));
    VDD site_pos_ = VDD(3, VD(2,0.0));
  
    // Binding probs as p_ij, (i for xlink, j for site)  
    VDD bind_probs_ = VDD(xlink_pos_.size(), VD(site_pos_.size(), 0.0));
    VDD curr_bind_probs_ = VDD(xlink_pos_.size(), VD(site_pos_.size(), 0.0));
    VDD unbind_probs_ = VDD(xlink_pos_.size(), VD(site_pos_.size(), 0.0));
    double tot_bind_prob_ = 0.0;
    
    // binding trackers
    VBB is_bound_ = VBB(xlink_pos_.size(), VB(site_pos_.size(), false));
    VBB bound_this_dt_ = VBB(xlink_pos_.size(), VB(site_pos_.size(), false));
    VI site_bonds_= VI(site_pos_.size(),0);
    VI xlink_bonds_ = VI(xlink_pos_.size(),0);
  
    // Binding trackers for determining collision
    VI site_bonds_prev_ = VI(site_pos_.size(),0);
    VI xlink_bonds_prev_ = VI(xlink_pos_.size(),0);

  public:
    Simulation();
    Simulation(double dt, int simType);
    ~Simulation();
    void Bind(int i, int j); // bind xlink i to site j
    void Unbind(int i, int j); // unbind xlink i from site j
    void Knockout(mt19937_64 &rng, uniform_real_distribution<double> &ran, bool i_wise);
    void Init();
    void Run(); // run once
    void Run(int n); // run n times
    template <class T>
    void Print(vector<vector<T> > &v);

    // Output vectors
    VDD bound_time = VDD(xlink_pos_.size(), VD(site_pos_.size(), 0.0));
    VDD bound_time_avg = VDD(xlink_pos_.size(), VD(site_pos_.size(), 0.0));
    VDD bound_time_95 = VDD(xlink_pos_.size(), VD(site_pos_.size(), 0.0));
    
    // 0 for collisions within a timestep allowed
    // 1 for first come, first serve within a timestep
    // 2 for random shuffling
    // 3 for whether to bind, then which
    // 4 for knockout method
    const int simType;
    
    // booleans derived from simType parameter
    bool allow_collisions_; // Allow collisions in main timestep
};

Simulation::Simulation(): dt_(1e-4), simType(0) { Init(); }
Simulation::Simulation(double dt, int type): dt_(dt), simType(type) { Init(); }
Simulation::~Simulation() { 
  ofile_.close();
}

template<class T>
void Simulation::Print(vector<vector<T> > &v) {
  iLOOP {
    jLOOP {
      cout << v[i][j] << " ";
      ofile_ << v[i][j];
      if ((i < 1) || (j < 2)) ofile_ << ",";
    }
    cout << endl;
  }
  ofile_ << endl;
}

// Perform knockout operation by comparing probabilities of each potential 
// conflict. i_wise decides whether to perform knockouts by first looping through
// xlinks (i_wise = true) or sites (i_wise = false)
void Simulation::Knockout(mt19937_64 &rng, uniform_real_distribution<double> &ran, 
                          bool i_wise) {
  int i, j;
  int *slow_ptr = &i; // points to slow/outer index (i if i_wise)
  int *fast_ptr = &j; // points to fast/inner index (j if i_wise)
  int slow_max = xlink_pos_.size();
  int fast_max = site_pos_.size();
  
  // swap slow, fast index if j_wise
  if (!i_wise) {
    slow_ptr = &j;
    fast_ptr = &i;
    slow_max = site_pos_.size();
    fast_max = xlink_pos_.size();
  }

  // Increment slow pointer index first, then fast pointer index
  for (*slow_ptr = 0; *slow_ptr < slow_max; (*slow_ptr)++) { 
    
    // Sum up conflicts
    double prob_conflicts = 0.0;
    int bound_count_ = 0;
    for (*fast_ptr = 0; *fast_ptr < fast_max; (*fast_ptr)++) {
      if (bound_this_dt_[i][j]) {
        prob_conflicts += bind_probs_[i][j];
        bound_count_++;
      }
    }
 
    // There is no conflict if just one is bound
    if (bound_count_ <= 1) continue;
    double roll_knockout = ran(rng);
    double cumul_prob = 0.0;
    int winner = -1;
    for (*fast_ptr = 0; *fast_ptr < fast_max; (*fast_ptr)++) {
      if (bound_this_dt_[i][j]) {
        cumul_prob += bind_probs_[i][j]/prob_conflicts;
        if (roll_knockout < cumul_prob) {
           winner = *fast_ptr; // Want inner index regardless of i_wise/j_wise
           break;
        }
      }
    }
    if (winner == -1) cerr << "knockout did not find winner!" << endl;
    
    //Unbind the loser sites/xlinks
    for (*fast_ptr = 0; *fast_ptr < fast_max; (*fast_ptr)++) {
      if (bound_this_dt_[i][j] && (*fast_ptr!=winner)) Unbind(i, j);
    }
  }
}

// Bind xlink i to site j
void Simulation::Bind(int i, int j) {
  is_bound_[i][j] = true;
  site_bonds_[j]++;
  xlink_bonds_[i]++;
  bound_this_dt_[i][j] = true;
}

// Unbind xlink i to site j
void Simulation::Unbind(int i, int j) {
  is_bound_[i][j] = false;
  site_bonds_[j]--;
  xlink_bonds_[i]--;
  bound_this_dt_[i][j] = false;
}

// Initialize simulation by filling vectors, calculating probs
void Simulation::Init() {
  // Open the output file
  string ofile_name = "output" + to_string(simType) + "_" + to_string(dt_) 
                      + ".csv";
  ofile_.open(ofile_name);

  // Set crosslink positions
  xlink_pos_[0][0] = 0.9; // x
  xlink_pos_[0][1] = 0.0; // y
  xlink_pos_[1][0] = 0.7; // x
  xlink_pos_[1][1] = 0.2; // y

  // Set site positions (10 degrees seperated on unit circle)
  jLOOP {
    site_pos_[j][0] = cos(10.0*j*M_PI/180.0); // x
    site_pos_[j][1] = -sin(10.0*j*M_PI/180.0); // y
  }
  
  // populate bind_probs_ and unbind_probs_. Ratio between bind prob and unbind prob must be
  // the Boltzmann factor.
  iLOOP jLOOP {
    bind_probs_[i][j] = kon_*dt_*exp(0.5*(E_0_-k_*(pow(xlink_pos_[i][0]-site_pos_[j][0],2)+
                        pow(xlink_pos_[i][1]-site_pos_[j][1],2))));
    unbind_probs_[i][j] = kon_*dt_*exp(-0.5*(E_0_-k_*(pow(xlink_pos_[i][0]-site_pos_[j][0],2)+
                          pow(xlink_pos_[i][1]-site_pos_[j][1],2))));
  }

  // derive booleans from simType
  // Allow collisions in Knockout type because we will knock them out afterwards
  allow_collisions_ = ((simType == 0) || (simType == 4));
}

// Run simulation once
void Simulation::Run() {
  
  // seed rng- using chrono and random libraries in order to get a [0,1) range and to
  // increase precision. Copy-pasted code from StackOverflow :)
  mt19937_64 rng;
  uint64_t timeSeed = chrono::high_resolution_clock::now().time_since_epoch().count();
  seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
  rng.seed(ss);
  uniform_real_distribution<double> ran(0.0, 1.0);
  int i_ind,j_ind;

  // Binding/unbinding main timestep loop
  for (double t = 0; t < (tmax_ - TOLERANCE); t += dt_) {

    // Initialize vectors for "whether/which" simTypes
    if (simType == 3) {
      iLOOP jLOOP curr_bind_probs_[i][j] = 0.0;
      tot_bind_prob_ = 0.0;
    }

    // Initialize bound-this-timestep vector for knockout simType
    if (simType == 4) {
      iLOOP jLOOP bound_this_dt_[i][j] = false;
    }

    // Shuffle ordering vectors for random shuffle simType
    if (simType == 2) {
      shuffle(xlink_index_.begin(), xlink_index_.end(), rng);
      shuffle(site_index_.begin(), site_index_.end(), rng);
    }

    iLOOP jLOOP {
      // Create random numbers and swap indices with shuffled vectors
      double roll_ij = ran(rng);
      i_ind = xlink_index_[i]; 
      j_ind = site_index_[j];
        
      // Bind if site and xlink were not bound in previous timestep (for coll's)
      // or if not currently bound (for all other simTypes)
      if ((allow_collisions_ && !site_bonds_prev_[j_ind] && !xlink_bonds_prev_[i_ind]) ||
          (!allow_collisions_ && !site_bonds_[j_ind] && !xlink_bonds_[i_ind])) {

        // If doing whether/which, save potential probabilities
        if (simType == 3) {
          tot_bind_prob_ += bind_probs_[i_ind][j_ind];
          curr_bind_probs_[i_ind][j_ind] = bind_probs_[i_ind][j_ind];
        
        // If not doing whether/which, directly roll to bind
        } else if (roll_ij < bind_probs_[i_ind][j_ind]) Bind(i_ind, j_ind);

      // Roll to unbind
      } else if (is_bound_[i_ind][j_ind] && (roll_ij < unbind_probs_[i_ind][j_ind])) {
        Unbind(i_ind, j_ind);
      }
    }

    // Make whether to bind/which to bind decision if using whether/which simType
    if (simType == 3) {
      double roll_whether_ = ran(rng);
      
      // Exponent for Poisson process
      if (roll_whether_ < (1-exp(-tot_bind_prob_))) {
        double roll_which = ran(rng);
        double cumul_prob = 0.0;
        [&] {
          iLOOP jLOOP {
            cumul_prob += curr_bind_probs_[i][j]/tot_bind_prob_;
            if (roll_which < cumul_prob) {
              Bind(i, j);
              return; // to end of iLOOP jLOOP
            }
            else if ((cumul_prob-1.0) > -TOLERANCE) {
              cerr << "error: which/whether run should never get here" << endl;
            }
          }
        }();
      }
    }

    if (simType == 4) {
      Knockout(rng, ran, true); // Knockout i_wise
      Knockout(rng, ran, false); // Knockout j_wise
    }

    // Discover if there are collisions, save new bound data
    iLOOP {
      if ((xlink_bonds_[i]-xlink_bonds_prev_[i])>1) {
        collision_count_++;
      }
      xlink_bonds_prev_[i] = xlink_bonds_[i];
    }
    jLOOP {
      if ((site_bonds_[j]-site_bonds_prev_[j])>1) {
        collision_count_++;
      }
      site_bonds_prev_[j] = site_bonds_[j];
    }
    // Record dt to bound percent after timestep
    iLOOP jLOOP if (is_bound_[i][j]) bound_time[i][j]+=dt_/tmax_;
  }
}

// Run n times
void Simulation::Run(int n) {
  VDDD bound_times = VDDD(n, VDD(xlink_pos_.size(), VD(site_pos_.size(), 0.0)));
  trialLOOP {
    Run();
    iLOOP jLOOP bound_time_avg[i][j] += bound_time[i][j]/n;
    bound_times[trial] = bound_time;

    // Set all trackers back to 0/false
    is_bound_ = VBB(xlink_pos_.size(), VB(site_pos_.size(), false));
    site_bonds_= VI(site_pos_.size(), 0);
    xlink_bonds_ = VI(xlink_pos_.size(), 0);
    site_bonds_prev_ = VI(site_pos_.size(), 0);
    xlink_bonds_prev_ = VI(xlink_pos_.size(), 0);
    bound_time = VDD(xlink_pos_.size(), VD(site_pos_.size(), 0.0));
  }
  trialLOOP iLOOP jLOOP {
    bound_time_95[i][j] += pow(bound_times[trial][i][j]-bound_time_avg[i][j], 2);
  }
  // 95% confidence interval is 2*sigma
  iLOOP jLOOP bound_time_95[i][j] = 2*pow(bound_time_95[i][j]*(1.0/n)*(1.0/(n-1)), 0.5);
    
  // Print out all info
  cout << "Collisions: " << collision_count_ << endl;
  cout << "Means:" << endl;
  Print(bound_time_avg);
  cout << endl << "+/-:" << endl;
  Print(bound_time_95);
}

int main() {
  double dt = 1e-4;

  // Run five trial simulations
  cout << endl << "Collisions allowed:" << endl;
  Simulation simAllowed(dt, 0);
  simAllowed.Run(10000);

  cout << endl << "FCFS:" << endl;
  Simulation simFCFS(dt, 1);
  simFCFS.Run(10000);

  cout << endl << "Sites/xlinks randomly shuffled:" << endl;
  Simulation simRand(dt, 2);
  simRand.Run(10000);

  cout << endl << "Whether/which:" << endl;
  Simulation simWW(dt, 3);
  simWW.Run(10000);

  cout << endl << "Knockout:" << endl;
  Simulation simKO(dt, 4);
  simKO.Run(10000);
  return 0;
}

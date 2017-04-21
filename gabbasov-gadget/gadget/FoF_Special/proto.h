

double periodic(double x);
double  periodic_wrap(double x);

void save_groups(char *indexlist_fname, char *catalogue_fname, char *particles_fname);
void save_particles(char *fname);
void reordering(void);
int course_binning(void);
void find_groups(void);
void check_cell(int p, int i, int j, int k);

int find_files(char *fname);
void loadpositions(char *fname, int files);

void allocate_memory(void);
void sort2_int(unsigned long n, int arr[], int brr[]);
void sort_int(unsigned long n, int arr[]);
 
void link_nearest(void);
void  determine_nearest(void);
int bin_for_nearsearch(void);

void sort2_fltint(unsigned long n, float arr[], int brr[], int crr[]);
void peano_hilbert_order(void);
int force_treebuild(void);
void force_update_node_recursive(int no, int sib, double len, double cx, double cy, double cz);
float ngb_select_closest(int k, int n, float *arr);
float ngb_treefind(float xyz[3], int desngb, int *pnearest, float hguess);

void reorder_particles(void);
void ngb_treeallocate(int maxnodes);
void ngb_treebuild(void);
int ngb_treefind_variable(float searchcenter[3], double hguess);



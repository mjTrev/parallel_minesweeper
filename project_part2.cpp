

#include <iostream>
#include <mpi.h>
#include <ctime>
#include <cstdlib>
#include <random>
#define MCW MPI_COMM_WORLD

using namespace std;

class Cell {
public:
	int row;
	int col;
	bool is_bomb=false;
	bool is_revealed=false;
	bool is_flagged=false;
	bool is_done = false;
	int neighbor_bombs = 0;
	int flagged_neighbors = 0;
	int revealed_neighbors = 0;
	int neighbor_count = 0;

	Cell(){
		row = 0;
		col = 0;
	}

	Cell(int i, int j){
		row = i;
		col = j;
	}
};

bool game_over = false;
int DIMENSIONS = 80;
int total_bombs = DIMENSIONS*DIMENSIONS*.2;
Cell neighbors[8];
Cell GRID[120][120];
int b_grid[120][120];
int f_grid[120][120];

void get_all_neighbors(int row, int col) {
	int i = 0;

	if (row - 1 >= 0){
		if (col - 1 >= 0) {
			neighbors[i] = GRID[row - 1][col - 1];
			i++;
		}
		if (col + 1 < DIMENSIONS) {
			neighbors[i] = GRID[row - 1][col + 1];
			i++;
		}
		neighbors[i] = GRID[row - 1][col];
		i++;
	}
	if (row + 1 < DIMENSIONS){
		if (col - 1 >= 0) {
			neighbors[i] = GRID[row + 1][col - 1];
			i++;
		}
		if (col + 1 < DIMENSIONS) {
			neighbors[i] = GRID[row + 1][col + 1];
			i++;
		}
		neighbors[i] = GRID[row + 1][col];
		i++;
	}
	if (col - 1 >= 0) {
		neighbors[i] = GRID[row][col - 1];
		i++;
	}
	if (col + 1 < DIMENSIONS) {
		neighbors[i] = GRID[row][col + 1];
		i++;
	}
	
	GRID[row][col].neighbor_count = i;
}

void fill_grid_data() {
	for (int row = 0; row < DIMENSIONS; row++) {
		for (int col = 0; col < DIMENSIONS; col++) {
			get_all_neighbors(row, col);
			for (int n = 0; n < GRID[row][col].neighbor_count; n++) {
				if (GRID[neighbors[n].row][neighbors[n].col].is_bomb) GRID[row][col].neighbor_bombs++;
			}
		}
	}
}

void display_solved_grid() {
	for (int i = 0; i < DIMENSIONS; i++) cout << "-";
	cout << endl;
	for (int row = 0; row < DIMENSIONS; row++) {
		for (int col = 0; col < DIMENSIONS; col++) {
			if (GRID[row][col].is_bomb) cout << "X";
			else cout << GRID[row][col].neighbor_bombs;
		}
		cout << endl;
	}
	for (int i = 0; i < DIMENSIONS; i++) cout << "-";
	cout << endl;
}

void display_grid() {
	for (int i = 0; i < DIMENSIONS; i++) cout << "-";
	cout << endl;
	for (int row = 0; row < DIMENSIONS; row++) {
		for (int col = 0; col < DIMENSIONS; col++) {
			if (GRID[row][col].is_revealed) {
				if (GRID[row][col].is_bomb) cout << "X";
				else cout << GRID[row][col].neighbor_bombs;
			}
			else if (GRID[row][col].is_flagged) {
				cout << "F";
			}
			else {
				cout << " ";
			}
		}
		cout << endl;
	}
	for (int i = 0; i < DIMENSIONS; i++) cout << "-";
	cout << endl;
}

void click(int row, int col) {
	
	if (GRID[row][col].is_revealed) return;
	GRID[row][col].is_revealed = true;
	
	get_all_neighbors(row, col);
	for (int i = 0; i < GRID[row][col].neighbor_count; i++) {
		GRID[neighbors[i].row][neighbors[i].col].revealed_neighbors++;
	}
	
	if (GRID[row][col].is_bomb) {
		game_over = true;
		return;
	}
	
	if (GRID[row][col].neighbor_bombs == 0) {
		for (int i = 0; i < GRID[row][col].neighbor_count; i++) {
			get_all_neighbors(row, col);
			click(neighbors[i].row, neighbors[i].col);
		}
	}
}

void flag(int row, int col) {
	if (GRID[row][col].is_flagged) return;
	GRID[row][col].is_flagged = true;

	get_all_neighbors(row, col);
	for (int i = 0; i < GRID[row][col].neighbor_count; i++) {
		GRID[neighbors[i].row][neighbors[i].col].flagged_neighbors++;
	}
}

bool make_move(int min_dim, int max_dim) {
	bool made_change = true;
	while (made_change) {
		made_change = false;
		for (int row = min_dim; row < max_dim; row++) {
			for (int col = min_dim; col < max_dim; col++) {
				//not done and not revealed
				if (!GRID[row][col].is_done && GRID[row][col].is_revealed) {
					//if no bomb neighbors or all neighbors are revealed or flagged then done
					if (GRID[row][col].neighbor_bombs == 0) GRID[row][col].is_done = true;
					else if (GRID[row][col].neighbor_count == GRID[row][col].revealed_neighbors + GRID[row][col].flagged_neighbors) GRID[row][col].is_done = true;
					else { //if not done
						//if unrevealed neighbors are all bombs
						if (GRID[row][col].neighbor_count - GRID[row][col].revealed_neighbors - GRID[row][col].is_flagged == GRID[row][col].neighbor_bombs) {
							for (int i = 0; i < GRID[row][col].neighbor_count; i++) {
								get_all_neighbors(row, col);
								if (!GRID[neighbors[i].row][neighbors[i].col].is_revealed) flag(neighbors[i].row, neighbors[i].col);
								made_change = true;
							}
						}
						if (GRID[row][col].flagged_neighbors == GRID[row][col].neighbor_bombs) {//if all bombs are flagged
							for (int i = 0; i < GRID[row][col].neighbor_count; i++) {
								get_all_neighbors(row, col);
								if (!GRID[neighbors[i].row][neighbors[i].col].is_flagged) click(neighbors[i].row, neighbors[i].col);
								made_change = true;
							}
						}
					}
				}
			}
		}
	}
}

void random_click(int min_dim, int dim) {
	srand(time(NULL));
	int row;
	int col;
	do {
		row = min_dim + rand() % dim;
		col = min_dim + rand() % dim;
	} while (GRID[row][col].is_revealed || GRID[row][col].is_flagged || GRID[row][col].is_bomb);

	click(row, col);
}

bool solved_section(int min_dim, int max_dim) {
	for (int row = min_dim; row < max_dim; row++) {
		for (int col = min_dim; col < max_dim; col++) {
			//all flags are on a bomb
			if (GRID[row][col].is_bomb && !GRID[row][col].is_flagged) return false;
			//all bombs are flagged
			if (!GRID[row][col].is_bomb && GRID[row][col].is_flagged) return false;
			//no bombs are revealed
			if (GRID[row][col].is_bomb && GRID[row][col].is_revealed) return false;
		}
	}
	game_over = true;
	return true;
}

void get_flags(int min_dim, int max_dim) {
	for (int i = min_dim; i < max_dim; i++) {
		for (int j = min_dim; j < max_dim; j++) {
			if (GRID[i][j].is_flagged) f_grid[i][j] = 1;
		}
	}
}

int main(int argc, char **argv){
	clock_t start;
	double duration;
	srand(time(NULL));
	int rand_num = 0;
	int rand_num2 = 0;
	int bomb_count = 0;
	int random_clicks = 0;

	//mpi stuff
	int rank, size;
	int data = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MCW, &rank);
	MPI_Comm_size(MCW, &size);
	int dim = DIMENSIONS / (size - 1);
	// initialize grid
	for (int i = 0; i < DIMENSIONS; i++) {
		for (int j = 0; j < DIMENSIONS; j++) {
			GRID[i][j] = Cell(i, j);
			b_grid[i][j] = 0;
			f_grid[i][j] = 0;
		}
	}

	//initialize 
	if (!rank) {
		cout << "start" << endl;
		start = clock();

		// fill some cells with bombs
		while (bomb_count < total_bombs) {
			rand_num = 0 + rand() % DIMENSIONS;
			rand_num2 = 0 + rand() % DIMENSIONS;

			if (!GRID[rand_num][rand_num2].is_bomb) {
				GRID[rand_num][rand_num2].is_bomb = true;
				b_grid[rand_num][rand_num2] = 1;
				bomb_count++;
			}
		}
		//send bomb list to other processes
		for (int i = 1; i < size; i++)
		{
			MPI_Send(&b_grid, 90000, MPI_INT, i, 0, MCW);
		}
		
		bool changes_made = true;
		//while(changes_made)
		//recieve code
		int code = 0;
		for (int i = 1; i < size; i++)
		{
			//all sub sections
			MPI_Recv(&data, 1, MPI_INT, MPI_ANY_SOURCE, 0, MCW, MPI_STATUSES_IGNORE);
			
		}



		duration = (clock() - start) / (double)CLOCKS_PER_SEC;
		cout << "Random clicks: " << random_clicks << endl;
		cout << duration << endl;
		
	}
	else{
		//put bombs in personal grid
		MPI_Recv(&b_grid, 90000, MPI_INT, 0, 0, MCW, MPI_STATUSES_IGNORE);
		for (int i = 0; i < DIMENSIONS; i++) {
			for (int j = 0; j < DIMENSIONS; j++) {
				if (b_grid[i][j]) GRID[i][j].is_bomb == true;
			}
		}
		
		// fill grid with data
		fill_grid_data();
		
		//start solving
		int min_dim = (rank - 1) * dim;
		int max_dim = min_dim + dim - 1;

		// first random click
		int first_row;
		int first_col;
		do {
			first_row = min_dim + dim / 3 + rand() % dim / 3;
			first_col = min_dim + dim / 3 + rand() % dim / 3;
		} while (GRID[first_row][first_col].is_bomb);
		click(first_row, first_col);

		
		//TODO: Simulate taking turns
		int status = 0;
		status = make_move(min_dim, max_dim);
		if (solved_section(min_dim, max_dim)) status = 3;
		MPI_Send(&status, 3, MPI_INT, 0, 0, MCW);
		if (status != 3) {
			//repeat
		}
	}
	
	MPI_Finalize();
	return 0;
}


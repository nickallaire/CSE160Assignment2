
double cput_time(void);
int usageRow();
int usageCol();
int verifyData(int argc, char c1[], char c2[], char c3[], char c4[], char c5[], char c6[], char c7[]);
double **calcArraySize(int row, int col);
double **fillInteriorCells(int startI, int startJ, int row, int col, double mean, double **array);
double **copyArray(int row, int col, double **arrayFrom, double **arrayTo);
double calcDiff(int row, int col, double **array, double **arrayOld, double diff);


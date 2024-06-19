#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <tuple>
#include <chrono>

#define E 2.718281828459045
#define PI 3.141592653589793





// global variables

static double dH = 0;
static double dT = 0;
static int N = 0;
static bool FILE_LOGGING = 0;

static double EPS = 1;

static double SIG = 0;

typedef struct Vector2
{
    double x = 0;
    double y = 0;
};

// Functions prototypes
//Exacts solution
double exata(double x, double y, double t);
double exataConv(double x, double y, double t);
Vector2 exataVec(double x,double y,double t );



//Debugging and exporting
void exportField(const std::vector<std::vector<double>> &field, const std::string &filename);
void exportVectorField(const std::vector<std::vector<Vector2>> &field, const std::string& filename);
void updateErrorCSV(const std::string &filename, int N, double dH, double dT, double error);



//Initialization
std::vector<std::vector<double>> initTriDiagonalMatrix(int size);
std::vector<double> initFontVector(std::vector<std::vector<double>> uAnt, int size, double initCondStart, double initCondEnd);



//Solutions
std::vector<std::vector<double>> SOLVE_ADI();
std::vector<std::vector<double>> SOLVE_ADI_CONV();
std::vector<std::vector<double>> SOLVE_EXACT();
std::vector<std::vector<double>> SOLVE_EXACT_CONV();
std::vector<std::vector<Vector2>> SOLVE_EXACT_VEC();

Vector2 VelocityFieldFunction(double x, double y);

//External algorithms
std::vector<double> TDMA(std::vector<double> lower, std::vector<double> main, std::vector<double> upper, std::vector<double> font);



//Utility

void printVector(std::vector<double> vet);
void printScalarField(std::vector<std::vector<double>> field);
void printVectorField(std::vector<std::vector<Vector2>> field);
std::vector<double> copy(std::vector<double> dest, std::vector<double> src, int offset);


const double TFINAL = 1.0;
const double TINITIAL = 0.0;



int main(int argc, char *argv[])
{

    N = std::stoi(argv[1]); // elementos na malha
    FILE_LOGGING = std::stoi(argv[2]);
    std::cout << "N = " << argv[1] << std::endl;
    dH = 1.0 / ((double)N - 1);
    dT = dH; //(dH); // Um passo no tempo com N = 4

    SOLVE_EXACT_VEC();
    /*

    SIG = (EPS * dT) / (dH * dH); // redefinido no convectivo  pois o eps muda

    std::cout << "dH = " << dH << std::endl;
    std::cout << "dT = " << dT << std::endl;
    std::cout << "Size of double: " << sizeof(double) << std::endl;
    // defining velocity field

    std::vector<std::vector<Vector2>> velocityField(N);

    // Exporting the vector field
    for (int i = 0; i < N; i++)
    {
        velocityField.push_back(std::vector<Vector2>());
        for (int j = 0; j < N; j++)
        {
            velocityField[i].push_back(VelocityFieldFunction(j * dH, i * dH));
        }
    }

    exportVectorField(velocityField,"VectorFields/VectorField_0.csv");

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<double>> uAprox = SOLVE_ADI_CONV();//SOLVE_ADI();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time = end - start;

    std::vector<std::vector<double>> uExact = SOLVE_EXACT_CONV();//SOLVE_EXACT();

    if (FILE_LOGGING)
    {
        std::cout << "Time for aprox solutions WITH exporting: " << time.count() << " seconds" << std::endl;
    }
    else
    {
        std::cout << "Time for aprox solutions WITHOUT exporting: " << time.count() << " seconds" << std::endl;
    }

    double errorMax = abs(uExact[0][0] - uAprox[0][0]);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (abs(uExact[i][j] - uAprox[i][j]) > errorMax)
            {
                errorMax = abs(uExact[i][j] - uAprox[i][j]);
            }
        }
    }

    std::cout << "Max error found: " << errorMax << std::endl;

    if (FILE_LOGGING)
        updateErrorCSV("Error/erros_dif_sqrt.csv", N, dH, dT, errorMax);

    std::cout << "Finished Simulation" << std::endl;

    */
}

std::vector<std::vector<double>> SOLVE_ADI()
{

    std::cout << "Started solving ADI..." << std::endl;
    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> uMid;
    std::vector<std::vector<double>> uAnt;
    int IT = 0;

    // Initial condition
    for (int i = 0; i < N; i += 1)
    {
        std::vector<double> u_row;
        std::vector<double> zero_row;

        for (int j = 0; j < N; j += 1)
        {

            u_row.push_back(exata(j * dH, i * dH, 0)); // cuidado com a ordem aqui, o eixo Y cresce pra baixo
            zero_row.push_back(0);
        }
        uMid.push_back(zero_row);
        u.push_back(zero_row);
        uAnt.push_back(u_row);
    }
    IT++;

    double time = 0;
    time += dT;
    exportField(uAnt, "ScalarFields/Field_0.csv");

    std::vector<std::vector<double>> matrix = initTriDiagonalMatrix(N - 2);

    // martrix é U[Y][X]

    while (time < TFINAL)
    {

        std::cout << "\r"
                  << "Solving IT " << IT << std::flush;

        //  Solving lines in X
        for (int i = 1; i < N - 1; i++)
        {
            // Imposting init conds
            uMid[i][0] = exata(0, i * dH, time);
            uMid[i][N - 1] = exata((N - 1) * dH, i * dH, time);

            int c = 0;
            std::vector<double> font(N - 2);
            for (int j = 1; j < N - 1; j++)
            {
                font[c] = 0.5 * SIG * (uAnt[i - 1][j] - (2 * uAnt[i][j]) + uAnt[i + 1][j]) + uAnt[i][j]; // add aqui a font
                c++;
            }

            font[0] += 0;               // exata(0, i * dH, time);                    // u[0][j];             // u[i][0];
            font[font.size() - 1] += 0; // exata(i * dH, i * dH, time); // u[i][N - 1];
            // std::cout << "Line " << i << " |";
            // printVector(font);

            std::vector<double> line = TDMA(matrix[0], matrix[1], matrix[2], font);

            c = 0;
            for (int k = 1; k < N - 1; k++)
            {
                uMid[i][k] = line[c];
                c++;
            }
        }

        time += (dT / 2);

        // Solving collums in Y
        for (int j = 1; j < N - 1; j++)
        {
            // Imposing init cond
            u[0][j] = exata(j * dH, 0, time); // The Y cresce pra baixo
            u[N - 1][j] = 0;

            std::vector<double> font(N - 2);
            int c = 0;
            for (int i = 1; i < N - 1; i++)
            {
                font[c] = 0.5 * SIG * (uMid[i][j + 1] - 2 * uMid[i][j] + uMid[i][j - 1]) + uMid[i][j];
                c++;
            }
            // font[0] += u[0][j];
            // font[font.size() - 1] += u[N - 1][j];
            font[0] += exata(j * dH, (N - 1) * dH, time);
            font[font.size() - 1] += exata(j * dH, 0, time);

            std::vector<double> collum = TDMA(matrix[0], matrix[1], matrix[2], font);

            // setting as colunas
            c = 0;
            for (int k = 1; k < N - 1; k++)
            {
                u[k][j] = collum[c];
                c++;
            }
        }

        // Save

        if (FILE_LOGGING)
        {
            std::stringstream fileName;
            fileName << "ScalarFields/Field_" << IT << ".csv";
            exportField(u, fileName.str());
        }
        IT++;

        // Reset buffers
        uAnt = u;

        time += (dT / 2);
    }
    std::cout << std::endl;

    std::cout << "Last Iteration stopped at: " << time << std::endl;
    return u;
}

std::vector<std::vector<double>> SOLVE_EXACT_CONV()
{
    std::cout << "Solving exact..." << std::endl;
    std::vector<std::vector<double>> u;
    int IT = 0;
    for (double n = 0; n < TFINAL; n += dT)
    {
        std::cout << "\r"
                  << "Exact at IT " << IT << std::flush;

        u.clear();
        for (double i = 0; i <= 1; i += dH)
        {
            std::vector<double> u_row;
            for (double j = 0; j <= 1; j += dH)
            {
                u_row.push_back(exataConv(j, i, n));
            }
            u.push_back(u_row);
        }

        if (FILE_LOGGING)
        {
            std::stringstream fileName;
            fileName << "ScalarFields/Field_" << IT << ".csv";
            exportField(u, fileName.str());

            fileName.clear();
        }
        IT++;
    };
    std::cout << std::endl;

    return u;
}

std::vector<std::vector<double>> SOLVE_EXACT()
{
    std::cout << "Solving exact..." << std::endl;
    std::vector<std::vector<double>> u;
    int IT = 0;
    for (double n = 0; n < TFINAL; n += dT)
    {
        std::cout << "\r"
                  << "Exact at IT " << IT << std::flush;

        u.clear();
        for (double i = 0; i <= 1; i += dH)
        {
            std::vector<double> u_row;
            for (double j = 0; j <= 1; j += dH)
            {
                u_row.push_back(exata(j, i, n));
            }
            u.push_back(u_row);
        }

        if (FILE_LOGGING)
        {
            std::stringstream fileName;
            fileName << "Fields/Field_" << IT << ".csv";
            exportField(u, fileName.str());

            fileName.clear();
        }
        IT++;
    };
    std::cout << std::endl;

    return u;
}


std::vector<std::vector<Vector2>> SOLVE_EXACT_VEC(){
std::cout << "Solving exact..." << std::endl;
    dH = 2.0 / ((double)N - 1);

    std::vector<std::vector<Vector2>> u;
    int IT = 0;
    for (double n = 0; n < TFINAL; n += dT)
    {
        std::cout << "\r"
                  << "Exact at IT " << IT << std::flush;

        u.clear();
        for (double i = -1; i <= 1; i += dH)
        {
            std::vector<Vector2> u_row;
            for (double j = -1; j <= 1; j += dH)
            {
                u_row.push_back(exataVec(j, i, n));
            }
            u.push_back(u_row);
        }

        if (FILE_LOGGING)
        {
            std::stringstream fileName;
            fileName << "VectorFields/VectorField_" << IT << ".csv";
            exportVectorField(u, fileName.str());
            fileName.clear();
        }
        IT++;
    };
    std::cout << std::endl;

    return u;


}
// Retorna uma lista de 3 elementos, contendo as 3 diagonais, list[0] a inferior, list[1] a main e list[2] a de cima
std::vector<std::vector<double>> initTriDiagonalMatrix(int size)
{
    std::vector<double> mainDiag(size);
    std::vector<double> bottomDiag(size - 1);
    std::vector<double> upperDiag(size - 1);

    for (int i = 0; i < size - 1; i++)
    {
        mainDiag[i] = 1 + SIG;
        bottomDiag[i] = -SIG * 0.5;
        upperDiag[i] = -SIG * 0.5;
    }
    mainDiag[size - 1] = 1 + SIG;

    std::vector<std::vector<double>> matrix;
    matrix.push_back(bottomDiag);
    matrix.push_back(mainDiag);
    matrix.push_back(upperDiag);

    return matrix;
}

// 0 no bool pra matri em X, 1 pra bool pra matrix em y
std::vector<std::vector<double>> initTriDiagonalMatrixConv(int size, std::vector<Vector2> fieldLine, bool XorY)
{
    std::vector<double> mainDiag(size);
    std::vector<double> bottomDiag(size - 1);
    std::vector<double> upperDiag(size - 1);
    if (XorY == false)
    {

        for (int i = 0; i < size - 1; i++)
        {
            double ro1 = (fieldLine[i + 1].x * dT) / dH;
            mainDiag[i] = 1 + SIG;
            bottomDiag[i] = (-ro1 / 4) + (-SIG * 0.5);
            upperDiag[i] = (ro1 / 4) + (-SIG * 0.5);
        }
        mainDiag[size - 1] = 1 + SIG;

        std::vector<std::vector<double>> matrix;
        matrix.push_back(bottomDiag);
        matrix.push_back(mainDiag);
        matrix.push_back(upperDiag);

        return matrix;
    }
    else
    {
        for (int i = 0; i < size - 1; i++)
        {
            double ro2 = (fieldLine[i + 1].y * dT) / dH;
            mainDiag[i] = 1 + SIG;
            bottomDiag[i] = (-ro2 / 4) + (-SIG * 0.5); // check for order in Y
            upperDiag[i] = (ro2 / 4) + (-SIG * 0.5);
        }
        mainDiag[size - 1] = 1 + SIG;

        std::vector<std::vector<double>> matrix;
        matrix.push_back(bottomDiag);
        matrix.push_back(mainDiag);
        matrix.push_back(upperDiag);

        return matrix;
    }
}

std::vector<double> TDMA(std::vector<double> lower, std::vector<double> mainDiagonal, std::vector<double> upper, std::vector<double> rhs)
{
    int n = mainDiagonal.size();
    std::vector<double> solution(n);

    // Forward elimination
    for (int i = 1; i < n; ++i)
    {
        double factor = lower[i - 1] / mainDiagonal[i - 1];
        mainDiagonal[i] -= factor * upper[i - 1];
        rhs[i] -= factor * rhs[i - 1];
    }

    // Backward substitution
    solution[n - 1] = rhs[n - 1] / mainDiagonal[n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        solution[i] = (rhs[i] - upper[i] * solution[i + 1]) / mainDiagonal[i];
    }

    return solution;
}

std::vector<std::vector<double>> SOLVE_ADI_CONV()
{
    EPS = 0.01;
    SIG = (EPS * dT) / (dH * dH);

    std::cout << "Started solving ADI..." << std::endl;
    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> uMid;
    std::vector<std::vector<double>> uAnt;

    std::vector<std::vector<Vector2>> velocityField(N);

    // Exporting the vector field
    for (int i = 0; i < N; i++)
    {
        velocityField.push_back(std::vector<Vector2>());
        for (int j = 0; j < N; j++)
        {
            velocityField[i].push_back(VelocityFieldFunction(j * dH, i * dH)); // The y axis grows donwards in the matrix because of reasons
        }
    }
    int IT = 0;

    // Initial condition
    for (int i = 0; i < N; i += 1)
    {
        std::vector<double> u_row;
        std::vector<double> zero_row;

        for (int j = 0; j < N; j += 1)
        {

            u_row.push_back(exataConv(j * dH, i * dH, 0)); // cuidado com a ordem aqui, o eixo Y cresce pra baixo
            zero_row.push_back(0);
        }
        uMid.push_back(zero_row);
        u.push_back(zero_row);
        uAnt.push_back(u_row);
    }
    IT++;

    double time = 0;
    time += dT / 2;
    exportField(uAnt, "Fields/Field_0.csv");

    while (time < TFINAL)
    {

        std::cout << "\r"
                  << "Solving IT " << IT << std::flush;

        for (int i = 0; i < N; i++)
        {

            // bounds conditions
            u[0][i] = exataConv(i * dH, 0, time + dT / 2);
            u[N - 1][i] = exataConv(i * dH, (N - 1) * dH, time + dT / 2);
            u[i][0] = exataConv(0, i * dH, time + dT / 2);
            u[i][N - 1] = exataConv((N - 1) * dH, i * dH, time + dT / 2);
        }

        for (int i = 1; i < N - 1; i++)
        {
            std::vector<std::vector<double>> matrix = initTriDiagonalMatrixConv(N - 2, velocityField[i], 0);
            // Bonds em x no passo intermed
            uMid[i][0] = exataConv(0, i * dH, time);
            uMid[i][N - 1] = exataConv((N - 1) * dH, i * dH, time);

            int c = 0;
            std::vector<double> font(N - 2);

            for (int j = 1; j < N - 1; j++)
            {
                double b2 = velocityField[i][j].y;
                double ro2 = (b2 * dT) / dH;

                font[c] = (1 - SIG) * uAnt[i][j] + (((SIG / 2) - (ro2 / 4)) * uAnt[i + 1][j]) + (((SIG / 2) + (ro2 / 4)) * uAnt[i - 1][j]);
                c++;
            }

            double b0 = velocityField[i][1].x;
            double bn = velocityField[i][N - 2].x;
            // mainDiag[i] = 1 + SIG;
            // bottomDiag[i] = (-ro1 / 4) + (-SIG * 0.5);
            // upperDiag[i] = (ro1 / 4) + (-SIG * 0.5);
            double ro0 = (b0 * dT) / dH;
            double ron = (bn * dT) / dH;
            font[0] += exataConv(0, i * dH, time) * ((ro0 / 4) + (SIG / 2));                         // u[0][j];             // u[i][0];
            font[font.size() - 1] += exataConv((N - 1) * dH, i * dH, time) * ((-ron / 4) + SIG / 2); // exata(i * dH, i * dH, time); // u[i][N - 1];
            // std::cout << "Line " << i << " |";
            // printVector(font);

            std::vector<double> line = TDMA(matrix[0], matrix[1], matrix[2], font);

            c = 0;
            for (int k = 1; k < N - 1; k++)
            {
                uMid[i][k] = line[c];
                c++;
            }
        }

        time += (dT / 2);

        // Solving collums in Y
        for (int j = 1; j < N - 1; j++)
        {
            // Imposing init cond
            u[0][j] = exataConv(j * dH, 0, time); // The Y axis grows downwards here
            u[N - 1][j] = exataConv(j * dH, (N - 1) * dH, time);

            // extracting the current collum of the velocity field
            std::vector<Vector2> velocityLine;

            for (int k = 0; k < N; k++)
            {
                velocityLine.push_back(velocityField[k][j]);
            }
            std::vector<std::vector<double>> matrix = initTriDiagonalMatrixConv(N - 2, velocityLine, 1);

            std::vector<double> font(N - 2);
            int c = 0;
            for (int i = 1; i < N - 1; i++)
            {
                double b1 = velocityField[i][j].x;
                double ro1 = (b1 * dT) / dH;
                font[c] = (1 - SIG) * uMid[i][j] + (((SIG / 2) - (ro1 / 4)) * uMid[i][j + 1]) + (((SIG / 2) + (ro1 / 4)) * uMid[i][j - 1]);
                c++;
            }
            // font[0] += u[0][j];
            // font[font.size() - 1] += u[N - 1][j];
            double b0 = velocityField[1][j].y;
            double bn = velocityField[N - 2][j].y;

            font[0] += exataConv(j * dH, 0, time) * ((b0 * dT) / (4 * dH) + SIG / 2);
            font[font.size() - 1] += exataConv(j * dH, (N - 1) * dH, time) * (-(bn * dT) / (4 * dH) + SIG / 2);

            std::vector<double> collum = TDMA(matrix[0], matrix[1], matrix[2], font);

            // setting the collum
            c = 0;
            for (int k = 1; k < N - 1; k++)
            {
                u[k][j] = collum[c];
                c++;
            }
        }

        // Saving this iteration field

        if (FILE_LOGGING)
        {
            std::stringstream fileName;
            fileName << "Fields/Field_" << IT << ".csv";
            exportField(u, fileName.str());
        }
        IT++;

        // Reset buffers
        uAnt = u;

        time += (dT / 2);
    }
    std::cout << std::endl;

    std::cout << "Last Iteration stopped at: " << time << std::endl;
    return u;
}

// Usada para a comparação do error e na cond de contorno
double exata(double x, double y, double t)
{

    double a = pow(E, -2 * PI * PI * t);
    return a * sin((PI * x)) * sin((PI * y));
}

// Exata com convecção
double exataConv(double x, double y, double t)
{
    double xline = (x * cos(4 * t)) + (y * sin(4 * t));
    double yline = (-x * sin(4 * t)) + (y * cos(4 * t));

    double gamma = 0.1;
    double eps = 0.01; // need check later if we change the define

    double xc = 0.2;
    double yc = 0.0;
    double coef = (2 * pow(gamma, 2)) / ((2 * pow(gamma, 2)) + (4 * eps * t));

    double xsqr = pow((xline - xc), 2);
    double ysqr = pow((yline - yc), 2);

    double expTerm = exp(-(xsqr + ysqr) / (2 * pow(gamma, 2) + (4 * eps * t)));

    return (expTerm * coef);
}


Vector2 exataVec(double x,double y,double t ){
    Vector2 vec;
    vec.x = -cos(PI*x)*sin(PI*y)*exp((-2*(PI*PI)*t)/20);
    vec.y = sin(PI*x)*cos(PI*y)*exp((-2*(PI*PI)*t)/20);

    return vec;

}

void exportField(const std::vector<std::vector<double>> &field, const std::string &filename)
{
    std::ofstream outfile(filename);

    size_t numRows = field.size();

    for (size_t i = 0; i < numRows; ++i)
    {
        bool isFirstValue = true; // for checking if it is the last collum

        for (const double value : field[i])
        {
            if (!isFirstValue)
            {
                outfile << ",";
            }
            outfile << value;
            isFirstValue = false;
        }

        // Dont put line break if last
        if (i < numRows - 1)
        {
            outfile << "\n"; // End of row
        }
    }

    outfile.close();
}

void exportVectorField(const std::vector<std::vector<Vector2>> &field, const std::string& filename){

    std::ofstream outfile(filename);

    size_t numRows = N;

    for (size_t i = 0; i < numRows; ++i)
    {
        bool isFirstValue = true; // for checking if it is the last collum

    //wtf is this
        for (const Vector2 value : field[i])
        {
            if (!isFirstValue)
            {
                outfile << ",";
            }
            outfile << value.x << ";"<<value.y;
            isFirstValue = false;
        }

        if (i < numRows - 1)
        {
            outfile << "\n"; // End of row
        }
    }

    outfile.close();

}

void printVector(std::vector<double> vet)
{

    for (std::vector<double>::iterator i = vet.begin(); i != vet.end(); i++)
    {

        std::cout << *(i) << ' ';
    }
    std::cout << std::endl;
}

void printScalarField(std::vector<std::vector<double>> field)
{
    for (auto i = field.begin(); i != field.end(); i++)
    {
        printVector(*(i));
    }
}
void printVectorField(std::vector<std::vector<Vector2>> field)
{
    for (auto i = field.begin(); i != field.end(); i++)
    {
        for (auto j = i->begin(); j != i->end(); j++)
        {
            std::cout << "(" << j->x << "," << j->y << ") ";
        }
    }
}

Vector2 VelocityFieldFunction(double x, double y)
{

    Vector2 v;
    v.x = -4.0 * y;
    v.y = 4.0 * x;
    return v;
}

// Essa função ta uma merda
void updateErrorCSV(const std::string &filename, int N, double dH, double dT, double error)
{
    std::ifstream file(filename);
    std::vector<int> NValues;
    std::vector<double> dHValues;
    std::vector<double> dTValues;
    std::vector<double> errorValues;

    if (file.is_open())
    {
        std::string line;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            int csvN;
            double csvdH, csvdT, csvError;

            if (iss >> csvN >> csvdH >> csvdT >> csvError)
            {
                NValues.push_back(csvN);
                dHValues.push_back(csvdH);
                dTValues.push_back(csvdT);
                errorValues.push_back(csvError);
            }
        }
        file.close();
    }

    auto it = std::find(NValues.begin(), NValues.end(), N);

    if (it != NValues.end())
    {
        size_t index = std::distance(NValues.begin(), it);
        dHValues[index] = dH;
        dTValues[index] = dT;
        errorValues[index] += error;
    }
    else
    {
        NValues.push_back(N);
        dHValues.push_back(dH);
        dTValues.push_back(dT);
        errorValues.push_back(error);
    }

    std::vector<std::tuple<int, double, double, double>> sortedData;
    for (size_t i = 0; i < NValues.size(); ++i)
    {
        sortedData.emplace_back(NValues[i], dHValues[i], dTValues[i], errorValues[i]);
    }

    // Sort the vector by N
    std::sort(sortedData.begin(), sortedData.end(), [](const auto &lhs, const auto &rhs)
              { return std::get<0>(lhs) < std::get<0>(rhs); });

    std::ofstream outFile(filename, std::ios::app);
    if (outFile.is_open())
    {
        for (const auto &entry : sortedData)
        {
            outFile << std::get<0>(entry) << ',' << std::get<1>(entry) << ',' << std::get<2>(entry) << ',' << std::get<3>(entry) << '\n';
        }
        outFile.close();
    }
    else
    {
        std::cerr << "Unable to open the output file." << std::endl;
    }
}
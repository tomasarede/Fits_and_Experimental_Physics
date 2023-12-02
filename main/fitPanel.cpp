#include <iostream>
#include <vector>
#include <random>
#include <TApplication.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include "TRandom3.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include <sstream>
#include <iomanip>



using namespace std;


void PlotGaussiansWithData(){
    // Read data from file
    ifstream file("/home/tomasarede/Desktop/3ano/LFEA/Arede_Codigos_Gamma/analysis/main/data1.txt");
    if (!file.is_open())
    {
        cout << "Failed to open data file." << endl;
        return;
    }


 // Vetores para armazenar os dados
    vector<double> x, y, dx, dy;
    
    // Ler os dados do arquivo
    double xValue, dxValue, yValue, dyValue;
    while (file >> xValue >> dxValue >> yValue >> dyValue) {
        x.push_back(xValue);
        dx.push_back(dxValue);
        y.push_back(yValue);
        dy.push_back(dyValue);
        
            
    }

    // Fechar o arquivo de dados
    file.close();
    

    // Create a TGraphErrors from the data
    int numPoints = x.size();
    TGraphErrors *graph = new TGraphErrors(numPoints, &x[0], &y[0], /*&dx[0]*/ nullptr, &dy[0]);

    graph->SetTitle("Variacao da posicao da fonte 22Na ao longo do eixo phi=0 graus");
    graph->GetXaxis()->SetTitle("Theta (graus)");
    graph->GetYaxis()->SetTitle("Coincidencias (cts)");

    graph->GetXaxis()->CenterTitle(true);
    graph->GetYaxis()->CenterTitle(true);

    // Create a canvas to display the plot
    TCanvas *canvas = new TCanvas("canvas", "Gaussians with Data", 950, 600);
    //canvas->SetLogy();

    // Draw the data points
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.7);
    graph->Draw("AP");

    // Create and draw the three Gaussian functions
    TF1 *gaussian1 = new TF1("gaussian1", "[0]*TMath::Gaus(x, [1], [2])", -60, 60);
    TF1 *gaussian2 = new TF1("gaussian2", "[0]*TMath::Gaus(x, [1], [2])", -60, 60);
    TF1 *gaussian3 = new TF1("gaussian3", "[0]*TMath::Gaus(x, [1], [2])", -60, 60);
    TF1 *gaussian4 = new TF1("gaussian4", "[0]*TMath::Gaus(x, [1], [2])", -60, 60);
    TF1 *gaussian5 = new TF1("gaussian5", "[0]*TMath::Gaus(x, [1], [2])", -60, 60);
    TF1 *gaussian6 = new TF1("gaussian6", "[0]*TMath::Gaus(x, [1], [2])", -60, 60);
    TF1 *gaussian7 = new TF1("gaussian7", "[0]*TMath::Gaus(x, [1], [2])", -60, 60);
    TF1 *gaussian8 = new TF1("gaussian8", "[0]*TMath::Gaus(x, [1], [2])", -60, 60);
    /*TF1 *gaussian9 = new TF1("gaussian9", "[0]*TMath::Gaus(x, [1], [2])", 0, 1023);
    TF1 *gaussian10 = new TF1("gaussian10", "[0]*TMath::Gaus(x, [1], [2])", 0, 1023);
    TF1 *gaussian11 = new TF1("gaussian11", "[0]*TMath::Gaus(x, [1], [2])", 0, 1023);
    TF1 *gaussian12 = new TF1("gaussian12", "[0]*TMath::Gaus(x, [1], [2])", 0, 1023);
    TF1 *gaussian13 = new TF1("gaussian12", "[0]*TMath::Gaus(x, [1], [2])", 0, 1023);
    TF1 *gaussian14 = new TF1("gaussian12", "[0]*TMath::Gaus(x, [1], [2])", 0, 1023);*/



    




    //Amp linear com ar
    gaussian1->SetParameters(6270, -31.65, 6.178);
    gaussian2->SetParameters(5706, -23.58, 6.094);
    gaussian3->SetParameters(5626,-15.79, 5.805);
    gaussian4->SetParameters(2300,-4.98,5.345);
    gaussian5->SetParameters(2361,8.952,6.153);
    gaussian6->SetParameters(2281,16.94,6.339);
    gaussian7->SetParameters(2363,24.35,6.201);
    gaussian8->SetParameters(2300,32.15,6.215);
    /*gaussian9->SetParameters(160.5,430.7,44.46);
    gaussian10->SetParameters(302.2,487.5,16.31);
    gaussian11->SetParameters(138.6,537.1,22.83);
    gaussian12->SetParameters(198.5,600.2,19.17);
    gaussian13->SetParameters(230.8,685.2,20.82);
    gaussian14->SetParameters(109.3,875.3,18.97);*/




    gaussian1->SetLineColor(kOrange + 1);
    gaussian2->SetLineColor(kViolet + 1);
    gaussian3->SetLineColor(kMagenta + 1);
    gaussian4->SetLineColor(kTeal + 1);
    gaussian5->SetLineColor(kAzure + 3);
    gaussian6->SetLineColor(kSpring + 3);
    gaussian7->SetLineColor(kPink + 3);
    gaussian8->SetLineColor(kCyan - 1);

    /*gaussian9->SetLineColor(kRed+1);
    gaussian10->SetLineColor(kRed+1);
    gaussian11->SetLineColor(kRed+1);
    gaussian12->SetLineColor(kRed+1);
    gaussian13->SetLineColor(kRed+1);
    gaussian14->SetLineColor(kRed+1);*/









    gaussian1->Draw("same");
    gaussian2->Draw("same");
    gaussian3->Draw("same");
    gaussian4->Draw("same");
    gaussian5->Draw("same");
    gaussian6->Draw("same");
    gaussian7->Draw("same");
    gaussian8->Draw("same");
    /*gaussian9->Draw("same");
    gaussian10->Draw("same");
    gaussian11->Draw("same");
    gaussian12->Draw("same");
    gaussian13->Draw("same");
    gaussian14->Draw("same");*/








    // Create and draw the sum of the Gaussians
    //TF1 *gaussianSum = new TF1("gaussianSum", "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5]) + [6]*TMath::Gaus(x, [7], [8])", 500, 650);

    // Amp Janela, primeiro ficheiro gonc_5_7.asc
    //gaussianSum->SetParameters(1.41771e+02, 4.29178e+02, 1.36281e+01, 1.04254e+03, 4.57667e+02, 7.76765e+00, 5.79581e+03, 4.79441e+02, 5.13494e+00);
    
    // Amp linear com ar
    //gaussianSum->SetParameters(6.41878e+01, 5.85332e+02, 2.98325e+01, 1.33457e+03, 5.92834e+02, 1.18222e+01, 2.79386e+03, 6.02739e+02, 6.32587e+00);
    
    // Amp Linear sem ar
    //gaussianSum->SetParameters(5.43254e+02, 6.97870e+02, 3.33783e+00, 3.93950e+03, 7.05066e+02, 2.11394e+00, 2.03306e+04, 7.09964e+02, 1.52146e+00);
    
    //gaussianSum->SetLineColor(kPink-4);

    //gaussianSum->Draw("same");

    // Create a legend and add entries for the data, individual Gaussians, and the sum
    TLegend *legend = new TLegend(0.9,0.6,0.7,0.9);
    legend->AddEntry(graph, "Dados", "lep");
    legend->AddEntry(gaussian1, "r = (-4.00 +/- 0.05) cm ", "l");
    legend->AddEntry(gaussian2, "r = (-3.00 +/- 0.05) cm ", "l");
    legend->AddEntry(gaussian3, "r = (-2.00 +/- 0.05) cm ", "l");
    legend->AddEntry(gaussian4, "r = (-1.00 +/- 0.05) cm ", "l");
    legend->AddEntry(gaussian5, "r = (1.00 +/- 0.05) cm ", "l");
    legend->AddEntry(gaussian6, "r = (2.00 +/- 0.05) cm ", "l");
    legend->AddEntry(gaussian7, "r = (3.00 +/- 0.05) cm ", "l");
    legend->AddEntry(gaussian8, "r = (4.00 +/- 0.05) cm ", "l");
    

    //legend->AddEntry(gaussian9, "Fits com Prob. < 5%", "l");







    //legend->AddEntry(gaussianSum, "Soma dos Picos", "l");
    legend->Draw();

    // Save the plot as a PDF file
    canvas->SaveAs("plot.pdf");
}


void Histogram_fit_gauss(){
// Criação do objeto TApplication
    TApplication app("app", nullptr, nullptr);


 
    // Abrir o arquivo de dados
    ifstream inputFile("/home/tomasarede/Arede_Codigos_Gamma/analysis/main/data.txt");
    if (!inputFile.is_open()) {
        cout << "Erro ao abrir o arquivo de dados." << endl;
    }

    // Ler os dados do arquivo e armazená-los em um vetor
    
    vector<double> data;
    double value;
    while (inputFile >> value) {
        //if (value !=0){
            data.push_back(value);
        //}
    }
    inputFile.close();
    
    

    
    // Configuração do gerador de números aleatórios
    /*random_device rd;
    mt19937 generator(rd());
    normal_distribution<double> distribution(0.0, 1.0);

    // Criar os pontos aleatórios
    vector<double> data;
    int numPoints = 1000;
    for (int i = 0; i < numPoints; ++i) {
        double value = distribution(generator);
        data.push_back(value);
    }
    */

    // Criar o histograma dos dados
    int numBins = 10;
    double minValue = *min_element(data.begin(), data.end());
    double maxValue = *max_element(data.begin(), data.end());

    TH1F *histogram = new TH1F("histogram", "Histogram", numBins, minValue, maxValue);
    histogram->SetTitle("Histogram");
    histogram->GetXaxis()->SetTitle("Value");
    histogram->GetYaxis()->SetTitle("Frequency");
    histogram->SetLineColor(kBlue);
    histogram->SetFillColor(kAzure - 9); // Preencher a área do histograma com a cor azul

    for (const auto& entry : data) {
        histogram->Fill(entry);
    }

    // Exibir o histograma e o fit na janela "c1"
    histogram->Draw();

    // Aguardar interação do usuário
    app.Run();

    // Limpar memória
    delete histogram;

}

void Fit_with_x_and_y_error(){


    TApplication app("app", nullptr, nullptr);

    // Abrir o arquivo de dados
    ifstream inputFile("/home/tomasarede/Desktop/3ano/LFEA/Arede_Codigos_Gamma/analysis/main/data1.txt");
    if (!inputFile.is_open()) {
        cout << "Erro ao abrir o arquivo data.txt!" << endl;
    }

    // Vetores para armazenar os dados
    vector<double> x, y, dx, dy;
    
    // Ler os dados do arquivo
    double xValue, dxValue, yValue, dyValue;
    while (inputFile >> xValue >> dxValue >> yValue >> dyValue) {
            x.push_back(xValue);
            dx.push_back(dxValue);
            y.push_back(yValue);
            dy.push_back(dyValue);
        
            
    }

    // Fechar o arquivo de dados
    inputFile.close();

    // Create a TGraphErrors from the data
    int numPoints = x.size();
    TGraphErrors *graph = new TGraphErrors(numPoints, &x[0], &y[0], /*&dx[0]*/  nullptr, &dy[0] /*nullptr*/);

    // Create a canvas to display the graph
    TCanvas *canvas = new TCanvas("canvas", "Data Fitting", 800, 600);
    graph->Draw("AP"); // "AP" option to display both markers and error bars
    graph->SetTitle("Variacao da fonte no eixo y");
    graph->GetXaxis()->SetTitle("Posicao y (cm)");
    graph->GetYaxis()->SetTitle("Taxa de Coincidencias");

    graph->GetXaxis()->CenterTitle(true);
    graph->GetYaxis()->CenterTitle(true);
    graph->SetMarkerSize(4.0);


    // Launch the ROOT Fit Panel
    canvas->Update(); // Update the canvas to display the graph
    canvas->cd();
    graph->FitPanel();
    
    // Wait for the Fit Panel to be closed
    app.Run();  
}


void Grafico_3D(/*const string& filename*/){

    TApplication app("app", nullptr, nullptr);

  // Personalizar o estilo do ROOT
  gStyle->SetPalette(1); // Configurar a paleta de cores como "calor"

  // Vetores para armazenar os dados
  
  vector<double> x, y, z;

  TRandom3 random; // Gerador de números aleatórios

  // Criar pontos aleatórios
  const int numPoints = 100;
  for (int i = 0; i < numPoints; ++i) {
    double t = random.Uniform(0, 10); // Tempo
    double xVal = t;
    double yVal = random.Gaus(0, 1); // Posição em y (gaussiana)
    double zVal = random.Gaus(0, 1); // Posição em z (gaussiana)

    x.push_back(xVal);
    y.push_back(yVal);
    z.push_back(zVal);
   }

  // Criar o gráfico 3D
  TGraph2D *graph = new TGraph2D(numPoints, &x[0], &y[0], &z[0]);

  // Criar um canvas para exibir o gráfico
  TCanvas *canvas = new TCanvas("canvas", "Data Fitting", 800, 600);
  graph->Draw("SURF1"); // "SURF1" option to display the graph as a continuous surface

  app.Run();

/*

  // Ler os dados do arquivo
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "Error opening file: " << filename << std::endl;
  }


  double x_val, ex_val, y_val, ey_val, z_val, ez_val;
  while (file >> x_val >> y_val >> z_val) {
    x.push_back(x_val);
    y.push_back(y_val);
    z.push_back(z_val);
  }

  file.close();
  */


}

void integrator(){
        
    // Gaussian 1 parameters
    double constant1 = 1.41771e+02;
    double mean1 = 4.29178e+02;
    double sigma1 = 1.36281e+01;

    // Gaussian 2 parameters
    double constant2 = 1.04254e+03;
    double mean2 = 4.57667e+02;
    double sigma2 = 7.76765e+00;

    // Gaussian 3 parameters
    double constant3 = 5.79581e+03;
    double mean3 = 4.79441e+02;
    double sigma3 = 5.13494e+00;
    

    /*
    // Gaussian 1 parameters
    double constant1 = 1.46865e+02;
    double mean1 = 4.27294e+02;
    double sigma1 = 1.10506e+01;

    // Gaussian 2 parameters
    double constant2 = 1.04464e+03;
    double mean2 = 4.56804e+02;
    double sigma2 = 6.52156e+00;

    // Gaussian 3 parameters
    double constant3 = 6.02019e+03;
    double mean3 = 4.79626e+02;
    double sigma3 = 5.03484e+00;
    */

    // Create TF1 objects for each Gaussian
    TF1 *gaussian1 = new TF1("gaussian1", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0, 1023);
    gaussian1->SetParameters(constant1, mean1, sigma1);

    TF1 *gaussian2 = new TF1("gaussian2", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0, 1023);
    gaussian2->SetParameters(constant2, mean2, sigma2);

    TF1 *gaussian3 = new TF1("gaussian3", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0, 1023);
    gaussian3->SetParameters(constant3, mean3, sigma3);

    // Perform integration for each Gaussian
    double integral1 = gaussian1->Integral(0, 441);
    double integral2 = gaussian2->Integral(441, 467.8);
    double integral3 = gaussian3->Integral(467.8, 1023);
    
    double sumIntegral = integral1 + integral2 + integral3;

    // Print the integration results
    std::cout << "Integral of Gaussian 1: " << integral1 << std::endl;
    std::cout << "Integral of Gaussian 2: " << integral2 << std::endl;
    std::cout << "Integral of Gaussian 3: " << integral3 << std::endl;
    std::cout << "Total Integral: " << sumIntegral << std::endl;
    std::cout << "Probability of Gaussian 1: " << integral1/sumIntegral*100 << std::endl;
    std::cout << "Probability of Gaussian 2: " << integral2/sumIntegral*100 << std::endl;
    std::cout << "Probability of Gaussian 3: " << integral3/sumIntegral*100 << std::endl;

    // Clean up
    delete gaussian1;
    delete gaussian2;
    delete gaussian3;

}
void print(vector<double> const &input)
{
    for (int i = 0; i < input.size(); i++) {
        cout << input.at(i) << ' ';
    }
}
 
void JoelhoERetroCalculator(int indice){
        // Abrir o arquivo de dados
   ifstream inputFile("/home/tomasarede/Arede_Codigos_Gamma/analysis/main/decay_eu.txt");
    if (!inputFile.is_open()) {
        cout << "Erro ao abrir o arquivo data.txt!" << endl;
    }

    // Vetor para armazenar os dados
    vector<double> x;

    // Ler os dados do arquivo
    string line;
    while (getline(inputFile, line)) {
        istringstream iss(line);
        double xValue;
        if (iss >> xValue) {
            x.push_back(xValue);
        }
    }

    // Fechar o arquivo de dados
    inputFile.close();


    int tamanho = x.size();
    vector<double> desvio1;
    vector<double> desvio2;
    vector<double> energy_values = {
        //41.4526919401083,
    	92.5039821599235,
        //127.817776361899,
        172.481682064352,
        //249.735584581077,
        //349.767441860465,
        //433.074227460975,
        561.299776999044,
        687.29531697993,
        //777.769990442816,
        856.7760433259,
        //957.285759796113,
        //1092.67919719656,
        //1395.48263778273
    };

    cout << "Energia Experimental: "<< energy_values[indice] << endl;
    for(int ii = 0; ii < tamanho; ii++){
        double Eg1=x[ii];
        double Eretro = Eg1/(1+2*Eg1/511);
        double Ecompton=Eg1-Eretro;
        desvio1.push_back(abs(energy_values[indice]-Eretro));
        desvio2.push_back(abs(energy_values[indice]-Ecompton));
    }
    
    double desvio1_min= *min_element(desvio1.begin(), desvio1.end());
    double desvio2_min= *min_element(desvio2.begin(), desvio2.end());
    //cout << desvio1_min<<endl;
    //cout << desvio2_min<<endl;

    auto iter_1 = find(desvio1.begin(), desvio1.end(), desvio1_min);
    int indice_1 = distance(desvio1.begin(), iter_1);

    auto iter_2 = find(desvio2.begin(), desvio2.end(), desvio2_min);
    int indice_2 = distance(desvio2.begin(), iter_2);
    int pos1= indice_1+1;
    int pos2= indice_2+1;

    double Eg_f=x[indice_1];
    double Eretro_f = Eg_f/(1+2*Eg_f/511);
    double Ecompton_f=Eg_f-Eretro_f;


    double Eg_f2=x[indice_2];
    double Eretro_f2 = Eg_f2/(1+2*Eg_f2/511);
    double Ecompton_f2=Eg_f2-Eretro_f2;


    
    cout << "Para a energia do pico com " << energy_values[indice] << " a retrodifusão mais próxima é do valor da posição " << pos1 << " que corresponde ao valor " << x[indice_1] <<"."<<endl;
    cout << "O nuclideo da posição "<< pos1 << ":"<<endl;
    cout << "Retrodifusão: " << Eretro_f << ";" << endl;
    cout << "Joelho de Compton: " << Ecompton_f << "." << endl;
    cout << endl;

    cout << "Para a energia do pico com " << energy_values[indice] << "o joelho de compton mais próximo é do valor da posição " << pos2 <<" que corresponde ao valor "<< x[indice_2]<<"."<<endl;
    cout << "O nuclideo da posição "<< pos2 << ":"<<endl;
    cout << "Retrodifusão: " << Eretro_f2 << ";" << endl;
    cout << "Joelho de Compton: " << Ecompton_f2 << "." << endl;

}

void table(){
    /*
    Facilita imenso fazer tabelas.
    Coloquem os valores e as incertezas em colunas diferentes no excel.
    Copiem e colem aqui.
    Depois corram o codigo.
    Copiem e colem o output no excel.
    Copiem e colem no vosso table generator (Ex:https://www.tablesgenerator.com/).
    Esta feito!
    */
    std::ifstream inputFile("/home/tomasarede/Arede_Codigos_Gamma/analysis/main/table.txt"); // Nome do arquivo de entrada

    if (!inputFile) {
        std::cout << "Erro ao abrir o arquivo de entrada." << std::endl;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);

        double coluna1, coluna2;
        char delimiter;

        if (iss >> coluna1 >> delimiter >> coluna2) {
            std::cout << "$" << coluna1 << " \\pm " << coluna2 << "$" << std::endl;
        }
    }

    inputFile.close();
}

void FitData()
{
    // Create a TApplication for the ROOT graphics
    TApplication app("app", nullptr, nullptr);

    // Read data from file
    std::ifstream file("/home/tomasarede/Desktop/3ano/LFEA/Arede_Codigos_Gamma/analysis/main/data1.txt");
    if (!file.is_open())
    {
        std::cout << "Failed to open data file." << std::endl;
        return;
    }

    // Arrays to store data
    std::vector<double> xValues;
    std::vector<double> yValues;
    std::vector<double> xErrors;
    std::vector<double> yErrors;


    double x, y, xErr, yErr;
    while (file >> x >> xErr >> y >> yErr)
    {
        xValues.push_back(x);
        yValues.push_back(y);
        xErrors.push_back(xErr);
        //yErrors.push_back( pow(abs(x),0.5)+1);
		yErrors.push_back( yErr );
    }
    file.close();

    // Create a TGraphErrors from the data
    int numPoints = xValues.size();
    TGraphErrors *graph = new TGraphErrors(numPoints, &xValues[0], &yValues[0], nullptr, &yErrors[0]);

    graph->SetTitle("R(#theta) (cts/s)");
    graph->GetXaxis()->SetTitle("#theta (degrees)");
    graph->GetYaxis()->SetTitle("R(#theta) (cts/s)");
	graph->GetXaxis()->SetLimits(-15,15);
    graph->GetXaxis()->CenterTitle(true);
    graph->GetYaxis()->CenterTitle(true);

    // Create a canvas to display the graph
    TCanvas *canvas = new TCanvas("canvas", "Data Fitting", 800, 600);
    graph->Draw("AP"); // "AP" option to display both markers and error bars

    // Define a custom TF1 function

    TF1 *fitFunc = new TF1("fitFunc", "2*[0]*[0]*acos(fabs(2*(-4)-[1]*x))/(2*[0])-fabs(2*(-4)-[1]*x)/2*sqrt(4*[0]*[0]-(2*(-4)-[1]*x)*(2*(-4)-[1]*x))", -4, 4);	//TF1 *fitFunc = new TF1("fitFunc", " ((1+TMath::Sign(1, 2*[1] - TMath::Abs(x-[2])))/2) * [0] * ( (TMath::Pi()/2) * (1 - TMath::Cos([1]*TMath::Pi()/180)) - TMath::ATan(TMath::Sqrt(2) *(TMath::Abs(x-[2])*TMath::Pi()/180/2) ) / TMath::Sqrt(TMath::Abs(TMath::Cos((x-[2])*TMath::Pi()/180) - TMath::Cos(2*[1]*TMath::Pi()/180) )) + TMath::ATan(TMath::Sqrt(2) *(TMath::Abs(x-[2])*TMath::Pi()/180/2)*TMath::Cos([1]*TMath::Pi()/180) ) / TMath::Sqrt(TMath::Abs(TMath::Cos((x-[2])*TMath::Pi()/180) - TMath::Cos(2*[1]*TMath::Pi()/180))) * TMath::Cos([1]*TMath::Pi()/180) ) + TMath::Abs([3]) ");
    //TF1 *fitFunc = new TF1("fitFunc", "[0] * ([1]^2 * TMath::ACos((x - [2]) / [1]) - TMath::Abs(x - [2]) * TMath::Sqrt([1]^2 - (x - [2])^2))");
    // TF1 *fitFunc = new TF1("fitFunc", "((1+TMath::Sign(1, TMath::Cos(x) - TMath::Cos(2*[1]))/2) * [0] * ((TMath::Pi()/2) * (1 - TMath::Cos([1])) - TMath::ATan(TMath::Abs(x) / TMath::Sqrt(2 * TMath::Cos(x) - 2 * TMath::Cos(2*[1])) + TMath::ATan(TMath::Cos([1]) * TMath::Abs(x) / TMath::Sqrt(2 * TMath::Cos(x) - 2 * TMath::Cos(2*[1])) * TMath::Cos([1])) + [3]");


    // Set initial parameter values and names
    fitFunc->SetParameters(625, 15.644, -1.175, 0);
    fitFunc->SetParNames("A", "beta", "theta0", "shift");


    // Fit the graph to the custom function
    graph->Fit(fitFunc, "R"); // "R" option for fit range using the graph's x-axis range
	graph->FitPanel();

    // Wait for the Fit Panel to be closed
    app.Run();
}
int main() {

    //Fazer o fit do histograma
    //Histogram_fit_gauss();

    //Fazer o fit conforme uma funcao
    Fit_with_x_and_y_error();

    //FitData();
    //Fazer o fit conforme uma funcao mas em 3D

    //string filename = "/home/tomasarede/Arede_Codigos_Gamma/analysis/main/data2.txt";
    //Grafico_3D(/*filename*/);
    
    //table();
    //PlotGaussiansWithData();
    //JoelhoERetroCalculator(3    );
    /*JoelhoERetroCalculator(1);
    JoelhoERetroCalculator(2);
    JoelhoERetroCalculator(3);
    JoelhoERetroCalculator(4);
    JoelhoERetroCalculator(5);
    JoelhoERetroCalculator(6);
    JoelhoERetroCalculator(7);
    JoelhoERetroCalculator(8);
    JoelhoERetroCalculator(9);
    JoelhoERetroCalculator(10);
    JoelhoERetroCalculator(11);*/
    

    return 0;
}







// 1215Square_HEX.cpp : 定义控制台应用程序的入口点。
//


#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <utility>
#include <ctime>
#include <numeric> 
#include <algorithm>
#include <numeric> 
#include <fstream>
#include <string>
using namespace std;



#define e 2.72;
//  variables representing the window size
float window_width = 350.0;
float window_height = 350.0;

vector<float> row(100, -1);
vector<vector<float>> grid(100, row); //格子数组
vector<pair<float, float>> res;//pds采样点位置记录


vector<int> indexRow(100, -1);
vector<vector<int>> indexMap(100, indexRow); //记录格子内点的序号


float r; //定义标准点距
int u = 1;//定义单位距离
float d = 0.0; //格子边长
int n = 0; //领域范围



enum Sampling
{
	PDS, SquarePixel, Hexagon
};

//sampling point set
void pointSet(int sampling, float radius) {

	grid.assign(100, row);
	indexMap.assign(100, indexRow);
	res.clear();


	int num = 0;//命中的数量
	int reject = 0;//连续拒绝次数
	bool accept = true;
	int index = 0;
	switch (sampling)
	{
	case 0:
		r = radius;
		d = sqrt(2)*r / 2;
		cout << "SquarePixel, radius=" << r ;
		for (float x = 0; x<window_width; x += r)
		{
			for (float y = 0; y < window_height; y += r) {

				res.push_back(pair<float, float>(x, y));

				int xGrid = x / d;
				int yGrid = y / d;
			}
		}
		break;
	case 1:
	{
		r = radius;
		d = sqrt(2)*r / 2;
		cout << "Hexagon, radius=" << r ;
		for (float x = 0.0; x<window_width; x += (r * 3 / 2)) {
			index++;
			float y = 0.0;
			if (index % 2 == 1)  y = sqrt(3)*r / 2;
			for (; y < window_height; y += (sqrt(3)*r)) {

				res.push_back(pair<float, float>(x, y));
				

				num++;
			}
		}
		index = 0;
		for (float x = r / 2; x<window_width; x += (r * 3 / 2))
		{
			index++;
			float y = 0.0;
			if (index % 2 == 0)  y = sqrt(3)*r / 2;
			for (; y < window_height; y += (sqrt(3)*r)) {

				res.push_back(pair<float, float>(x, y));
				
				num++;
			}
		}
	}
	break;

	case 2:
		r = radius;
		d = sqrt(2)*r / 2;
		cout << "PDS, radius=" << r ;
		srand(static_cast <unsigned> (time(0)));

		while (reject <= 150000) {

			accept = true;
			float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_width);
			float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_height);

			int xGrid = x / d;
			int yGrid = y / d;

			for (int i = xGrid - 2; i <= xGrid + 2; i++) {
				if (i < 0) continue;
				for (int j = yGrid - 2; j <= yGrid + 2; j++)
				{
					if (j<0) continue;
					if (indexMap[i][j] != -1) {
						int index = indexMap[i][j];
						float dX = res[index].first - x;
						float dY = res[index].second - y;

						float dis = sqrt(pow(dX, 2) + pow(dY, 2));
						if (dis < r) {
							accept = false;
							break;
						}
					}
				}
				if (!accept) {
					reject++;
					break;
				}
			}
			if (accept)
			{
				reject = 0;
				indexMap[xGrid][yGrid] = num;

				num++;
				res.push_back(pair<float, float>(x, y));
			}
		}
		break;
	default:
		break;
	}
	
	cout << " success: " << num << endl;

}


int main()
{
	while (true)
	{
		int sampling;//采样方法，grid领域范围，循环次数
		float radius;
		vector<float> resMean;

		cout << "samoling method: ";
		cin >> sampling;
		if (sampling < 0 || sampling>2) {
			cout << "input invalid" << endl;
			cout << endl;
			continue;
		}

		cout << "radius: ";
		cin >> radius;
		if (radius < 8 || radius>13) {
			cout << "input invalid" << endl;
			cout << endl;
			continue;
		}

		
		pointSet(sampling, radius);
		
		string file = "sampling";
		file.append(to_string(res.size()));
		file.append(".txt");

		ofstream outfile;
		outfile.open(file);
		outfile << res.size() << endl;
		for (auto i :res)
		{
			outfile << i.first / window_height << " " << i.second / window_width << endl;
		}
		outfile.close();
	}
}









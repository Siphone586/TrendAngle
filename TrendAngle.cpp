// TrendAngle.cpp : 定义 DLL 应用程序的导出函数。
//
#define _CRT_SECURE_NO_WARNINGS

#include "stdafx.h"
#include "stdio.h"
#include <string.h>
//
#include"TrendAngle.h"
#include <algorithm>
#include <windows.h>
#include <math.h>
#define M_PI 3.1416
//输入 真实波幅，计算波幅变化趋势角度（核滤波）

unsigned LENGTH = 8;;
	unsigned SCALE=2;
	float SMOOTH=2.;
	float FACTOR=1.;



//1.计算角度，弧度
float atan2my(float y, float x) {
	float angle = 0.0;
	if (x > 0) {
		angle = atan(y / x);
	}
	else {
		if (x < 0 && y >= 0) {
			angle = atan(y / x) + M_PI;
		}
		else {
			if (x < 0 && y < 0) {
				angle = atan(y / x) - M_PI;
			}
			else {
				if (x == 0 && y > 0) {
					angle = M_PI / 2;
				}
				else {
					if (x == 0 && y < 0) {
						angle = -M_PI / 2;
					}
				}
			}
		}
	}
	return angle;
}

//2 弧度转角度
float rad2deg(float angle_rad) {
	return angle_rad * 180 / M_PI;
}
//3 Epanechnikov 核函数
float epanechnikov_kernel(float *src, int length, float smooth, float factor) {
	float current_weight = 0.0;
	float cumulative_weight = 0.0;
	for (int i =length-1; i>0; i--) {
		float y = src[i];
		float u = pow(i, 2) / (pow(smooth, 2) * factor);
		float w = (u >= 1) ? 0 : (3. / 4) * (1 - pow(u, 2));
		current_weight += y * w;
		cumulative_weight += w;
	}
	return cumulative_weight > 0 ? current_weight / cumulative_weight : 0; // To avoid division by zero.  
}


void AtrAngle1(int DataLen, float*pfOUT, float*pfINa, float*pfINb, float*pfINc)
{
	//Na :atr   Nb close  Nc Lenth*1000+scale, Na要拆分为 2个量 
	//_size: 核的长度。 _h : 核的宽度。 _r : 核的锐度。
	unsigned param = int(pfINa[DataLen - 1]);//组合参数，要拆分;
	LENGTH = int(param / 1000); //核长度
	param %= 1000; // 去掉已经提取的千分位数字  
	SCALE = param ;	// 提取百分位数字  
	// LENGTH = pfINa[0];
	//SCALE = pfINb[0];
	SMOOTH =pfINb[DataLen - 1];
	FACTOR=pfINc[DataLen - 1];
}	

//TrendAngle.cpp:定义DLL应用程序的导出函数。
//
// 2号函数，计算趋势角度
void AtrAngle2(int DataLen, float*pfOUT, float*pfINa, float*pfINb, float*pfINc)
{
	//Na :atr   Nb close  Nc NULL
	//_size: 核的长度。 _h : 核的宽度。 _r : 核的锐度。
	float *atr = pfINa;//波幅
	float *close= pfINb; //收盘价

	float *slope = new float[DataLen+1];
	float *angle_rad= new float[DataLen + 1];
	float *angle_deg = new float[DataLen + 1];
	float *degrees = new float[DataLen + 1];

	if (DataLen > (LENGTH + 1))
	{


		//float *atr = malloc(sizeof(float) * DataLen);

		  // 初始化变量
		for (int i = 0; i < DataLen; i++) {
			slope[i] = 0.0;
			angle_deg[i] = 0.0;
			degrees[i] = 0.0;
		}
		//大循环
		for (int i = LENGTH; i < DataLen; i++)
		{
			slope[i] = (close[i] - close[i - LENGTH]) / (atr[i] / (LENGTH / SCALE) * LENGTH);
			angle_rad[i] = atan2(slope[i], 1.0);
			angle_deg[i] = rad2deg(angle_rad[i]);
			if (i > LENGTH * 2)
			{
				degrees[i] = epanechnikov_kernel(&angle_deg[i-LENGTH], LENGTH, SMOOTH, FACTOR);
				pfOUT[i] = angle_deg[i];
			}
		}
	}
	// 释放内存
	delete[]angle_rad;
	delete[]slope;
	delete[]angle_deg;
	delete[]degrees;

	return ;

		



	//delete[]buf;
}

// 3号函数，计算趋势角度 核后
void AtrAngle3(int DataLen, float*pfOUT, float*pfINa, float*pfINb, float*pfINc)
{
	//Na :atr   Nb close  Nc NULL 
	//_size: 核的长度。 _h : 核的宽度。 _r : 核的锐度。
	float *atr = pfINa;//波幅
	float *close = pfINb; //收盘价


	float *slope = new float[DataLen + 1];
	float *angle_rad = new float[DataLen + 1];
	float *angle_deg = new float[DataLen + 1];
	float *degrees = new float[DataLen + 1];

	if (DataLen > (LENGTH + 1))
	{


		//float *atr = malloc(sizeof(float) * DataLen);

		  // 初始化变量
		for (int i = 0; i < DataLen; i++) {
			slope[i] = 0.0;
			angle_deg[i] = 0.0;
			degrees[i] = 0.0;
		}
		//大循环
		for (int i = LENGTH; i < DataLen; i++)
		{
			slope[i] = (close[i] - close[i - LENGTH]) / (atr[i] / (LENGTH / SCALE) * LENGTH);
			angle_rad[i] = atan2(slope[i], 1.0);
			angle_deg[i] = rad2deg(angle_rad[i]);
			if (i > LENGTH * 2)
			{
				degrees[i] = epanechnikov_kernel(&angle_deg[i - LENGTH], LENGTH, SMOOTH, FACTOR);
				pfOUT[i] = degrees[i];
			}
		}
	}
	// 释放内存
	delete[]angle_rad;
	delete[]slope;
	delete[]angle_deg;
	delete[]degrees;

	return;





	//delete[]buf;
}


// 4号函数，核滤波
void Epanech4(int DataLen, float*pfOUT, float*pfINa, float*pfINb, float*pfINc)
{
	//Na :atr   Nb close  Nc Lenth*1000+scale*100+smooth*10+factor, NC要拆分为 3个量 
	//_size: 核的长度。 _h : 核的宽度。 _r : 核的锐度。
	float *atr = pfINa;//波幅
	float *close = pfINb; //收盘价
	float *slope = new float[DataLen + 1];
	float *angle_rad = new float[DataLen + 1];
	float *angle_deg = new float[DataLen + 1];
	float *degrees = new float[DataLen + 1];

	if (DataLen > (LENGTH + 1))
	{


		//float *atr = malloc(sizeof(float) * DataLen);

		  // 初始化变量
		for (int i = 0; i < DataLen; i++) {
			slope[i] = 0.0;
			angle_deg[i] = 0.0;
			degrees[i] = 0.0;
		}
		//大循环
		for (int i = LENGTH; i < DataLen; i++)
		{
			//slope[i] = (atr[i] - atr[i - LENGTH]) / (LENGTH);
			//angle_rad[i] = atan2(slope[i], 1.0);
			//angle_deg[i] = rad2deg(angle_rad[i]);
			if (i > LENGTH * 2)
			{
				degrees[i] = epanechnikov_kernel(&atr[i - LENGTH], LENGTH, SMOOTH, FACTOR);
				pfOUT[i] = degrees[i];
			}
		}
	}
	// 释放内存
	delete[]angle_rad;
	delete[]slope;
	delete[]angle_deg;
	delete[]degrees;

	return;





	//delete[]buf;
}






//加载的函数

PluginTCalcFuncInfo g_CalcFuncSets[] =
{


	{ 2,(pPluginFUNC)&AtrAngle2 },//计算趋势角度    核后
	{ 1,(pPluginFUNC)&AtrAngle1 },//计算趋势角度
	{ 0,NULL },
};



//导出给TCalc的注册函数

BOOL RegisterTdxFunc(PluginTCalcFuncInfo**pFun)
{
	if (*pFun == NULL)
	{
		(*pFun) = g_CalcFuncSets;

		return TRUE;
	}
	return FALSE;
}



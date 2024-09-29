// TrendAngle.cpp : ���� DLL Ӧ�ó���ĵ���������
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
//���� ��ʵ���������㲨���仯���ƽǶȣ����˲���

unsigned LENGTH = 8;;
	unsigned SCALE=2;
	float SMOOTH=2.;
	float FACTOR=1.;



//1.����Ƕȣ�����
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

//2 ����ת�Ƕ�
float rad2deg(float angle_rad) {
	return angle_rad * 180 / M_PI;
}
//3 Epanechnikov �˺���
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
	//Na :atr   Nb close  Nc Lenth*1000+scale, NaҪ���Ϊ 2���� 
	//_size: �˵ĳ��ȡ� _h : �˵Ŀ�ȡ� _r : �˵���ȡ�
	unsigned param = int(pfINa[DataLen - 1]);//��ϲ�����Ҫ���;
	LENGTH = int(param / 1000); //�˳���
	param %= 1000; // ȥ���Ѿ���ȡ��ǧ��λ����  
	SCALE = param ;	// ��ȡ�ٷ�λ����  
	// LENGTH = pfINa[0];
	//SCALE = pfINb[0];
	SMOOTH =pfINb[DataLen - 1];
	FACTOR=pfINc[DataLen - 1];
}	

//TrendAngle.cpp:����DLLӦ�ó���ĵ���������
//
// 2�ź������������ƽǶ�
void AtrAngle2(int DataLen, float*pfOUT, float*pfINa, float*pfINb, float*pfINc)
{
	//Na :atr   Nb close  Nc NULL
	//_size: �˵ĳ��ȡ� _h : �˵Ŀ�ȡ� _r : �˵���ȡ�
	float *atr = pfINa;//����
	float *close= pfINb; //���̼�

	float *slope = new float[DataLen+1];
	float *angle_rad= new float[DataLen + 1];
	float *angle_deg = new float[DataLen + 1];
	float *degrees = new float[DataLen + 1];

	if (DataLen > (LENGTH + 1))
	{


		//float *atr = malloc(sizeof(float) * DataLen);

		  // ��ʼ������
		for (int i = 0; i < DataLen; i++) {
			slope[i] = 0.0;
			angle_deg[i] = 0.0;
			degrees[i] = 0.0;
		}
		//��ѭ��
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
	// �ͷ��ڴ�
	delete[]angle_rad;
	delete[]slope;
	delete[]angle_deg;
	delete[]degrees;

	return ;

		



	//delete[]buf;
}

// 3�ź������������ƽǶ� �˺�
void AtrAngle3(int DataLen, float*pfOUT, float*pfINa, float*pfINb, float*pfINc)
{
	//Na :atr   Nb close  Nc NULL 
	//_size: �˵ĳ��ȡ� _h : �˵Ŀ�ȡ� _r : �˵���ȡ�
	float *atr = pfINa;//����
	float *close = pfINb; //���̼�


	float *slope = new float[DataLen + 1];
	float *angle_rad = new float[DataLen + 1];
	float *angle_deg = new float[DataLen + 1];
	float *degrees = new float[DataLen + 1];

	if (DataLen > (LENGTH + 1))
	{


		//float *atr = malloc(sizeof(float) * DataLen);

		  // ��ʼ������
		for (int i = 0; i < DataLen; i++) {
			slope[i] = 0.0;
			angle_deg[i] = 0.0;
			degrees[i] = 0.0;
		}
		//��ѭ��
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
	// �ͷ��ڴ�
	delete[]angle_rad;
	delete[]slope;
	delete[]angle_deg;
	delete[]degrees;

	return;





	//delete[]buf;
}


// 4�ź��������˲�
void Epanech4(int DataLen, float*pfOUT, float*pfINa, float*pfINb, float*pfINc)
{
	//Na :atr   Nb close  Nc Lenth*1000+scale*100+smooth*10+factor, NCҪ���Ϊ 3���� 
	//_size: �˵ĳ��ȡ� _h : �˵Ŀ�ȡ� _r : �˵���ȡ�
	float *atr = pfINa;//����
	float *close = pfINb; //���̼�
	float *slope = new float[DataLen + 1];
	float *angle_rad = new float[DataLen + 1];
	float *angle_deg = new float[DataLen + 1];
	float *degrees = new float[DataLen + 1];

	if (DataLen > (LENGTH + 1))
	{


		//float *atr = malloc(sizeof(float) * DataLen);

		  // ��ʼ������
		for (int i = 0; i < DataLen; i++) {
			slope[i] = 0.0;
			angle_deg[i] = 0.0;
			degrees[i] = 0.0;
		}
		//��ѭ��
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
	// �ͷ��ڴ�
	delete[]angle_rad;
	delete[]slope;
	delete[]angle_deg;
	delete[]degrees;

	return;





	//delete[]buf;
}






//���صĺ���

PluginTCalcFuncInfo g_CalcFuncSets[] =
{


	{ 2,(pPluginFUNC)&AtrAngle2 },//�������ƽǶ�    �˺�
	{ 1,(pPluginFUNC)&AtrAngle1 },//�������ƽǶ�
	{ 0,NULL },
};



//������TCalc��ע�ắ��

BOOL RegisterTdxFunc(PluginTCalcFuncInfo**pFun)
{
	if (*pFun == NULL)
	{
		(*pFun) = g_CalcFuncSets;

		return TRUE;
	}
	return FALSE;
}



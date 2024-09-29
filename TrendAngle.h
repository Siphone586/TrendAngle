#ifndef __PLUGIN_TCALC_FUNC
#define __PLUGIN_TCALC_FUNC


#pragma pack(push,1) 


#pragma pack(push,1)

//函数(数据个数,输出,输入a,输入b,输入c)
typedef void(*pPluginFUNC)(int, float*, float*, float*, float*);

typedef struct tagPluginTCalcFuncInfo
{
	unsigned short                nFuncMark;//函数编号
	pPluginFUNC                        pCallFunc;//函数地址
}PluginTCalcFuncInfo;

typedef BOOL(*pRegisterPluginFUNC)(PluginTCalcFuncInfo**);

//RegisterTdxFunc函数为注册函数，用户的函数全部都放入PluginTCalcFuncInfo结构的全局数组g_CalcFuncSets中。
typedef BOOL(*pRegisterPluginFUNC)(PluginTCalcFuncInfo**);

//=======================================

extern"C"
{
	//__cplusplus
	__declspec(dllexport)BOOL RegisterTdxFunc(PluginTCalcFuncInfo**pFun);
}


#pragma pack(pop)
#endif
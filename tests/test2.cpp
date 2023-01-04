/*
 * @Author: your name
 * @Date: 2021-11-11 19:21:29
 * @LastEditTime: 2021-11-11 20:19:26
 * @LastEditors: Please set LastEditors
 * @Description: 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * @FilePath: /Solver/test2.cpp
 */
#include "commonType.h"
#include "Zonotope.h"
#include <iostream>

using namespace std;

int main(){
    reachSolver::Zonotope<double> myzono(10);
    //NonlinearSys<double> mysys(6);
    //NonlinearSys<double> mysys(FunctionTest, 6, 1, 6);
    //NonlinearSys<double> mysys();
}

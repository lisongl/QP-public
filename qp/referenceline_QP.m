function [referenceline_x,referenceline_y] = ...
    fcn(referenceline_x_init,referenceline_y_init,w_cost_smooth,w_cost_length,w_cost_ref,x_lb,x_ub,y_lb,y_ub)
%该函数将平滑referenceline
%使用二次规划
%0.5x'Hx + f'x = min
% lb < x < ub
% 输入 referenceline x,y,init未平滑的参考线
% w_cost_smooth,w_cost_length,w_cost_ref 平滑代价，紧凑代价，几何相似代价权重系数
% x_lb,x_ub,y_lb,y_ub 允许x,y变化的上下界
% 输出 referenceline_x,referenceline_y 平滑后的参考线

%在这里为了简化，不使用参考线拼接算法
%计算量问题通过调度解决，规划周期为100ms，只要在100ms计算完就可以
%如果到最后计算速度成问题在优化，目前暂时不使用拼接

%simulink无法直接使用自带的二次规划函数，需要从外部引入

coder.extrinsic("quadprog");
% 二次规划形式
% H1 = w_cost_smooth*(A1'*A1) + w_cost_length*(A2'*A2) + w_cost_ref*(A3'*A3)
% f = -2 * w_cost_ref * referenceline_init
% x'H1x + f'x = 0.5 * x'(2H1)*x + f'x
% A1 = [1  0 -2  0  1  0
%       0  1  0 -2  0  1
%             1  0 -2  0  1  0
%             0  1  0 -2  0  1
%                   ...............
%A2 = [1  0 -1  0
%      0  1  0 -1
%            1  0 -1  0
%            0  1  0 -1
%                  ...........
%A3 为单位矩阵
%设置求解矩阵规模，因为参考线初始点一共181个，每个点有两个坐标，所以求解矩阵规模是362个
n = 181;
%referenceline初始化
referenceline_x = zeros(n,1);
referenceline_y = zeros(n,1);

A1 = zeros(2 * n - 4, 2 * n);
A2 = zeros(2 * n - 2, 2 * n);
A3 = eye(2 * n, 2 * n);
f = zeros(2 * n, 1);
lb = zeros(2 * n, 1);
ub = zeros(2 * n, 1);
persistent is_first_run;
persistent pre_result_x;
persistent pre_result_y;
persistent pre_referenceline_x_init;
persistent pre_referenceline_y_init;

if isempty(is_first_run)
    is_first_run = 0;
    for i = 1 : n
        f(2 * i - 1) = referenceline_x_init(i);
        f(2 * i) = referenceline_y_init(i);
        lb(2 * i - 1) = f(2 * i - 1) + x_lb;
        ub(2 * i - 1) = f(2 * i - 1) + x_ub;
        lb(2 * i) = f(2 * i) + y_lb;
        ub(2 * i) = f(2 * i) + y_ub;
    end

    for j = 1 :2: 2 * n - 5
        A1(j , j) = 1;
        A1(j , j + 2) = -2;
        A1(j , j + 4) = 1;
        A1(j + 1, j + 1) = 1;
        A1(j + 1, j + 3) = -2;
        A1(j + 1, j + 5) = 1;
    end

    for k = 1 : 2: 2 * n - 3
        A2(k, k) = 1;
        A2(k, k + 2) = -1;
        A2(k + 1, k + 1) = 1;
        A2(k + 1, k + 3) = -1;
    end

    H1 = w_cost_smooth * (A1'*A1) + w_cost_length * (A2'*A2) + w_cost_ref * A3;
    H = 2 * H1;
    %用 referenceline 的初值作为二次规划的迭代起点
    X0 = f;
    f = -2 * w_cost_ref * f;
    %二次规划求解
    X = quadprog(H,f,[],[],[],[],lb,ub,X0);
    %将结果输出到referenceline中
    for i = 1 : n
        referenceline_x(i) = X(2 * i - 1);
        referenceline_y(i) = X(2 * i);
    end
    pre_result_x = referenceline_x;
    pre_result_y = referenceline_y;
    pre_referenceline_x_init = referenceline_x_init;
    pre_referenceline_y_init = referenceline_y_init;
else
    if pre_referenceline_x_init(1) == referenceline_x_init(1) && ...
            pre_referenceline_y_init(1) == referenceline_y_init(1)
        %起点相同，可以认为本周期和上周期的referenceline_init 一模一样，就直接复用上个周期的结果
        referenceline_x = pre_result_x;
        referenceline_y = pre_result_y;
    else
        %起点不同。重复第一次运行的逻辑
        for i = 1 : n
            f(2 * i - 1) = referenceline_x_init(i);
            f(2 * i) = referenceline_y_init(i);
            lb(2 * i - 1) = f(2 * i - 1) + x_lb;
            ub(2 * i - 1) = f(2 * i - 1) + x_ub;
            lb(2 * i) = f(2 * i) + y_lb;
            ub(2 * i) = f(2 * i) + y_ub;
        end

        for j = 1 : 2: 2 * n - 5
            A1(i , i) = 1;
            A1(i , i + 2) = -2;
            A1(i , i + 4) = 1;
            A1(i + 1, i + 1) = 1;
            A1(i + 1, i + 3) = -2;
            A1(i + 1, i + 5) = 1;
        end

        for k = 1 : 2: 2 * n - 3
            A2(i, i) = 1;
            A2(i, i + 2) = -1;
            A2(i + 1, i + 1) = 1;
            A2(i + 1, i + 3) = -1;
        end

        H1 = w_cost_smooth * (A1'*A1) + w_cost_length * (A2'*A2) + w_cost_ref * A3;
        H = 2 * H1;
        %用 referenceline 的初值作为二次规划的迭代起点
        X0 = f;
        f = -2 * w_cost_ref * f;
        %二次规划求解
        X = quadprog(H,f,[],[],[],[],lb,ub,X0);
        %将结果输出到referenceline中
        for i = 1 : n
            referenceline_x(i) = X(2 * i - 1);
            referenceline_y(i) = X(2 * i);
        end
        pre_result_x = referenceline_x;
        pre_result_y = referenceline_y;
        pre_referenceline_x_init = referenceline_x_init;
        pre_referenceline_y_init = referenceline_y_init;
        
    end
    
end











    
 











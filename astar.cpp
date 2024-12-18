
#include<global_planner/astar.h>
#include<costmap_2d/cost_values.h>

namespace global_planner {

AStarExpansion::AStarExpansion(PotentialCalculator* p_calc, int xs, int ys) :
        Expander(p_calc, xs, ys) {
}

bool AStarExpansion::calculatePotentials(unsigned char* costs, double start_x, double start_y, double end_x, double end_y,
                                        int cycles, float* potential) {
    queue_.clear();
    queue2_.clear();

    int start_i = toIndex(start_x, start_y);//把起点一维化
    queue_.push_back(Index(start_i, 0));//将起点加入队列1末尾
    int end_i = toIndex(end_x, end_y);//把终点一维化
    queue2_.push_back(Index(end_i, 1e10));//将终点加入队列2末尾

    std::fill(potential, potential + ns_, POT_HIGH);//初始化潜力图，全为1.0e10
    potential[start_i] = 0;//将起点设为0
    potential[end_i] = 1e10;//将终点设为1e10（极大，方便）

    int goal1_i = toIndex(end_x, end_y);
    int goal2_i = toIndex(end_x, end_y);
    int cycle = 0;

    while (queue_.size() > 0 && cycle < cycles) {
        Index top1 = queue_[0];
        std::pop_heap(queue_.begin(), queue_.end(), greater1());//优先队列
        queue_.pop_back();
        Index top2 = queue2_[0];
        std::pop_heap(queue2_.begin(), queue2_.end(), less1());//优先队列
        queue2_.pop_back();

        //起点开始的搜索搜到终点
        int i = top1.i;
        if (i == goal1_i)
            return true;
        
        //终点开始的搜索搜到起点
        int j = top2.i;
        if (j == goal2.i)
        {
            //清空queue_队列
            while(!queue_.empty()){
                queue_.pop_back();
            }
            while(!queue2_empty()){
                //取堆顶元素
                Index out = queue2_[0];
                std::pop_heap(queue2_.begin(), queue2_.end(), greater1());
                queue2_.pop_back();
                //倒序入堆
                queue_.push_back(out);
                std::push_heap(queue_.begin(), queue_.end(), less1());
            }
            return true;
        }

        //双向搜索相遇(next_i相同时)
        float before_cost = top1.cost; //起点开始的已有代价
        int i_next[4]={i+1,i-1,i+nx_,i-nx_};//点i的四周
        int j_next[4]={j+1,j-1,j+nx_,j-nx_};//点j的四周

        for(int i_dir=0;i_dir<4;i_dir++)
            for(int j_dir=0;j_dir<4;j_dir++){
                if(i_next[i_dir] == j_next[j_dir]){
                    keep_queue_(before_cost,i_next[i_dir]);
                    return true;
                }
            }

        //正向搜索
        add(costs, potential, potential[i], i + 1, end_x, end_y);
        add(costs, potential, potential[i], i - 1, end_x, end_y);
        add(costs, potential, potential[i], i + nx_, end_x, end_y);
        add(costs, potential, potential[i], i - nx_, end_x, end_y);
        //反向搜索
        add1(costs, potential, potential[j], j + 1, start_x, start_y);
        add1(costs, potential, potential[j], j - 1, start_x, start_y);
        add1(costs, potential, potential[j], j + nx_, start_x, start_y);
        add1(costs, potential, potential[j], j - nx_, start_x, start_y);

        cycle++;
    }

    return false;
}
//起点开始添加点
void AStarExpansion::add(unsigned char* costs, float* potential, float prev_potential, int next_i, int end_x,
                         int end_y) {
    if (next_i < 0 || next_i >= ns_) //边界点，不选择
        return;

    if (potential[next_i] < POT_HIGH) //障碍物，不选择
        return;

    if(costs[next_i]>=lethal_cost_ && !(unknown_ && costs[next_i]==costmap_2d::NO_INFORMATION)) 
    //代价大于致死成本且该点不是代价地图上无信息的未知的点 ，不选择
        return;

    //不是以上的情况，则下一节点相邻节点中预估代价最小节点的代价加上当前代价作为下一点的代价
    potential[next_i] = p_calc_->calculatePotential(potential, costs[next_i] + neutral_cost_, next_i, prev_potential);
    /*potential[next_i] = p_calc_->calculatePotential1(potential, - costs[next_i] - neutral_cost_, next_i, prev_potential);*/
    int x = next_i % nx_, y = next_i / nx_; //x与y坐标
    float distance = abs(end_x - x) + abs(end_y - y); //曼哈顿估计，仅计算水平与垂直方向方格的数量总和

    queue_.push_back(Index(next_i, potential[next_i] + distance * neutral_cost_));//压栈，压入下一路径点与其估计总代价（已行进代价与预估代价）
    std::push_heap(queue_.begin(), queue_.end(), greater1());//压栈后的优先队列从大到小排序，即当前的正向路径
}
//终点开始添加点
void AStarExpansion::add1(unsigned char* costs, float* potential, float prev_potential, int next_i, int end_x,
                         int end_y) {
    if (next_i < 0 || next_i >= ns_) //边界点，不选择
        return;

    if (potential[next_i] < POT_HIGH) //障碍物，不选择
        return;

    if(costs[next_i]>=lethal_cost_ && !(unknown_ && costs[next_i]==costmap_2d::NO_INFORMATION)) 
    //代价大于致死成本且该点不是代价地图上无信息的未知的点 ，不选择
        return;

    //不是以上的情况，则下一节点相邻节点中预估代价最小节点的代价加上当前代价作为下一点的代价
    potential[next_i] = p_calc_->calculatePotential1(potential, costs[next_i] + neutral_cost_, next_i, prev_potential);
    int x = next_i % nx_, y = next_i / nx_; //x与y坐标
    float distance = abs(end_x - x) + abs(end_y - y); //曼哈顿估计，仅计算水平与垂直方向方格的数量总和

    queue2_.push_back(Index(next_i, potential[next_i] - distance * neutral_cost_));//压栈，压入下一路径点与其估计总代价（已行进代价与预估代价）
    std::push_heap(queue2_.begin(), queue2_.end(), less1());//压栈后的优先队列从小到大排序，即当前的反向路径
}

void AStarExpansion::keep_queue_(float before_cost,int i){//维护queue_队列的完整性
            int cost_add = 1; //之后每个节点的代价都默认+1
            queue_.push_back(Index(i,before.cost + 1));
            
            while(!queue2_empty()){
                cost_add++;
                //取堆顶元素
                Index out = queue2_[0];
                std::pop_heap(queue2_.begin(), queue2_.end(), greater1());
                queue2_.pop_back();
                //入堆
                queue_.push_back(Index(out.i ,before.cost + cost_add));
                std::push_heap(queue_.begin(), queue_.end(), greater1());
            }
}

} //end namespace global_planner

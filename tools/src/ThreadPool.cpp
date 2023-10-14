/**
 * @file ThreadPool.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief  
 * @version udef
 * @date 2023/09/20
 * 
 * @copyright Copyright (c) 2023 XLAB
 * 
 * modification history :
 * Date:      Version:    Author:
 * Changes:
 */

#include "tools/ThreadPool.h"

namespace NSPthread {

using namespace std;

ThreadPool::ThreadPool(size_t threadCount) : terminate(false) {
    for(int i = 0; i<threadCount; i++) {
        threads.emplace_back(&ThreadPool::threadFunction, this);
        threadCount++;
    }
}

void ThreadPool::addTask(shared_ptr<Task> taskptr) {
    {
        lock_guard<mutex> lock(queueMutex);  // 等待互斥量解锁
        taskQueue.emplace_back(taskptr);  // 添加任务
        taskCount++;
    }
    condition.notify_one();  // 通知线程
}

void ThreadPool::threadFunction() {
    while(true) {
        shared_ptr<Task> pTask;  //抽象类 Task 无法实例化，只能用指针
        {
            unique_lock<mutex> lock(queueMutex);
            condition.wait(lock, [this]() {return !taskQueue.empty() || terminate;});

            if(terminate && taskQueue.empty()) {
                break;
            }

            pTask = taskQueue.front();
            taskQueue.pop_front();
        }
        pTask->run();  // execute the task.
        taskCount--;
        completedTaskCount++;
    }
}

ThreadPool::~ThreadPool() {
    terminate = true;
    condition.notify_all();  //唤醒所有等待中的线程

    for(thread& th : threads) {
        if (th.joinable()) {
            th.join(); // 等待线程执行完毕
        }
    } 
}
}
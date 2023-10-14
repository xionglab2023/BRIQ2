/**
 * @file ThreadPool.h
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

#ifndef TOOLS_THREADPOOL_H_
#define TOOLS_THREADPOOL_H_


#include <thread>
#include <vector>
#include <deque>
#include <condition_variable>
#include <functional>
#include <iostream>

namespace NSPthread {
    using namespace std;

    class Task {  // Task 基类
        private:
            int priority;
            static atomic<size_t> requestCount;
            atomic<bool> isCancelRequired;
        protected:
            size_t id;
            clock_t createTime;
        public:
            Task():isCancelRequired(false),createTime(clock()),id(requestCount++){};
            virtual ~Task(){};

            virtual int run() = 0;  //纯虚函数，需要子类实现

            virtual int onCanceled(){
                return 1;
            }

            virtual int onCompleted(int){
                return 1;
            }

            size_t getID() {
                return id;
            }

            bool getCancelRequired() {
                return isCancelRequired;
            }

            void setCancelRequired() {
                isCancelRequired = true;
            }

            int getPriority() const {
                return priority;
            }

    };

    atomic<size_t> Task::requestCount = 0;

    class IntFuncTask:public Task  // 执行一个具有整数返回值的函数的 Task
    {
        private:
            function<int(void)> wf;
        public:
            IntFuncTask(function<int(void)> f): wf(f) {}
            IntFuncTask(): wf(nullptr){}

            virtual ~IntFuncTask(){}

            template <typename F, typename... FArgs, typename... Args>
            void asynBind(F(*f)(FArgs...), Args... args)
            {
                wf = bind(f,args...);
            }

            virtual int run() throw()
            {
                if(wf == nullptr)
                    return 86;
                try {
                    return wf();
                } catch(exception& ex) {
                    cerr << ex.what() <<endl;
                    return EXIT_FAILURE;
                } catch(const string& estr) {
                    cerr << estr << endl;
                    return EXIT_FAILURE;
                }
            }
    };

    class ThreadPool {
        private:
            void threadFunction();  //线程执行函数
            vector<thread> threads;
            // 任务队列相关
            deque<shared_ptr<Task> > taskQueue;             // 任务队列
            mutex queueMutex;                  // 任务队列访问互斥量
            condition_variable condition;      // 任务队列条件变量，通知线程有新任务可执行

            atomic<bool> terminate;  // 线程池终止标记

            //统计信息
            atomic<size_t> threadCount;  // 线程数量
            atomic<size_t> taskCount;  // 任务数量
            atomic<size_t> completedTaskCount;  //已完成任务数量
            chrono::steady_clock::time_point startTime;  // 线程池启动时间
 
        public:
            ThreadPool(size_t nt);
            void addTask(shared_ptr<Task> taskptr);  // 任务队列管理
            void addThreads(size_t nt);
            virtual ~ThreadPool();

            size_t getThreadCount() const {
                return threadCount.load();
            }

            size_t getTaskCount() const {
                return taskCount.load();
            }

            size_t getCompletedTaskCount() const {
                return completedTaskCount.load();
            }

            double getRunningTimeInSeconds() const {
                chrono::duration<double> duration = chrono::steady_clock::now() - startTime;
                return duration.count();
            }
    };
}

#endif /* TOOLS_THREADPOOL_H_ */
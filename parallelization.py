import multiprocessing


class Parallelization:

    @staticmethod
    def parallelize_7(func, *args_list):
        if len(args_list[0][0]) > 7:
            rest = len(args_list[0][0]) % 7
            last = 0
            for i in range(0,len(args_list[0][0]),7):
                if i + 6 < len(args_list[0][0]):
                    jobs = []
                    print('parallelization in 7 cores')
                    for cpu in range(7):
                        last = i+cpu
                        print('cpu '+str(cpu))
                        args = []
                        for param in args_list[0]:
                            args.append(param[i+cpu])
                        p = multiprocessing.Process(target=func(*args))
                        jobs.append(p)
                        p.start()
                else:
                    remaining = (last+rest)-(last)
                    print('rest of parameters: '+str(remaining)+' - parallelization in '+str(remaining)+' cores')
                    jobs = []
                    for cpu in range(remaining):
                        print('cpu '+str(cpu))
                        args = []
                        for param in args_list[0]:
                            args.append(param[i+cpu])
                        p = multiprocessing.Process(target=func(*args))
                        jobs.append(p)
                        p.start()
        else:
            print('size of parameters: '+str(len(args_list[0][0]))+' - parallelization in '+str(len(args_list[0][0]))+' cores')
            jobs = []
            for cpu in range(len(args_list[0][0])):
                print('cpu '+str(cpu))
                args = []
                for param in args_list[0]:
                    args.append(param[cpu])
                p = multiprocessing.Process(target=func(*args))
                jobs.append(p)
                p.start()
import multiprocessing
import db_handling
from progress.bar import Bar
import signal
import sys

class Parallelization:

    stop_process = False

    def signal_handler(self, signal, frame):
        self.stop_process = True

    def parallelize_7(self, func, *args_list, log_param=False, bar=False, handle_signal=False):
        if handle_signal:
            signal.signal(signal.SIGINT, self.signal_handler)

        protDB = db_handling.ProteinDatabase()

        if len(args_list[0][0]) > 7:
            rest = len(args_list[0][0]) % 7
            last = 0
            for i in range(0,len(args_list[0][0]),7):
                if i + 6 < len(args_list[0][0]):
                    jobs = []
                    # print('parallelization in 7 cores')
                    for cpu in range(7):
                        last = i+cpu
                        # print('cpu '+str(cpu))
                        args = []
                        for param in args_list[0]:
                            args.append(param[i+cpu])
                            # print(args)
                        p = multiprocessing.Process(target=func(*args))
                        jobs.append(p)
                        p.start()
                        if bar:
                            bar.next()
                else:
                    remaining = (last+rest)-(last)
                    # print('rest of parameters: '+str(remaining)+' - parallelization in '+str(remaining)+' cores')
                    jobs = []
                    for cpu in range(remaining):
                        # print('cpu '+str(cpu))
                        args = []
                        for param in args_list[0]:
                            args.append(param[i+cpu])
                            # print(args)
                        p = multiprocessing.Process(target=func(*args))
                        jobs.append(p)
                        p.start()
                        if bar:
                            bar.next()
                if log_param:
                    protDB.update_state({'comparisons_latest_seq_a':args[0][0], 'comparisons_latest_seq_b':args[1][0]},log_param)
                if handle_signal and self.stop_process:
                    sys.exit('\nprocess terminated correctly')
        else:
            # print('size of parameters: '+str(len(args_list[0][0]))+' - parallelization in '+str(len(args_list[0][0]))+' cores')
            jobs = []
            for cpu in range(len(args_list[0][0])):
                # print('cpu '+str(cpu))
                args = []
                for param in args_list[0]:
                    args.append(param[cpu])
                p = multiprocessing.Process(target=func(*args))
                jobs.append(p)
                p.start()
                if bar:
                    bar.next()
                if log_param:
                    protDB.update_state({'comparisons_latest_seq_a':args[0][0], 'comparisons_latest_seq_b':args[1][0]},log_param)
                if self.stop_process:
                    sys.exit('process terminated correctly')
                


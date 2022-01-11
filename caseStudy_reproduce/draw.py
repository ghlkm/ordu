import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

def readdata(filename):
    f = open(filename, "r")
    data = []
    str = f.readline()
    while str:
        str = str.strip()
        ar = str.split(' ')
        ar = [float(i.strip()) for i in ar]
        data.append([(ar[1] + ar[3]) / 2, (ar[2] + ar[4]) / 2])
        str = f.readline()
    f.close()
    return data

def readoption(filename):
    f = open(filename, "r")
    data = dict()
    str = f.readline()
    while str:
        str = str.strip()
        ar = str.split(',')
        data[ar[0].strip()]=[int(i.strip()) for i in ar[1:]]
        str = f.readline()
    f.close()
    return data

def main():
    datafile = str(sys.argv[1])
    optionfile = str(sys.argv[2])

    data=readdata(datafile)
    option_id=readoption(optionfile)

    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    plt.rc('font', family='serif')

    fig = plt.figure(figsize=(10, 10))
    ax = fig.gca()

    data = np.array(data)
    option=dict()
    s=set()
    for key in option_id:
        option_id[key]=[i-1 for i in option_id[key]]
        option[key]=np.array([data[i] for i in option_id[key]])
        for id in option_id[key]:
            s.add(id)
    rest = np.array([i for i in range(len(data)) if i not in s])
    ax.plot(data[rest, 0], data[rest, 1], '.', c='gray', markersize=10)

    # support maximum 4 kinds of results
    for key, mk, cl in zip(option, [6, 4, 5, 7], ['r', 'g', 'b', 'orange']):
        ax.plot(option[key][:, 0], option[key][:, 1], marker=mk, c=cl, markersize=15, linestyle='None', label=key)
    ax.legend(prop={'size': 28}, loc=4, ncol=2, columnspacing=0, handletextpad=0, borderpad=0.15)
    ax.set_xlabel('Points', fontsize=28)
    ax.set_ylabel('Rebounds', fontsize=28)
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)

    plt.show()
    # plt.savefig('.\\k2m6_pts_trb_43_57-4.pdf',format='pdf', bbox_inches='tight')



if __name__=='__main__':
    main()


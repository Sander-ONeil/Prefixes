import numpy as np


def vec(x):
    return np.array(x,dtype =np.float64)
def inv(x):
    return np.linalg.inv(x)


mlta = 5


prefixes ={
    'Mol':'6.02214076e+23',
    'T': '1000000000000',
    'G': '1000000000',
    'M': '1000000',
    'k': '1000',
    '':1,
    'c': '0.01',
    'm': '0.001',
    'μ': '0.000001',
    'n': '0.000000001',
    'p': '0.000000000001',
    'y': '0.000000000000000000000001'
    
}
for p in range(len(prefixes)):
    ind = list(prefixes)[p]
    prefixes[ind] = float(prefixes[ind])


#Base units are M L T A temp(K)

units =  {
    'g':vec([1,0,0,0,0]),
    'N':vec([1,1,-2,0,0]),
    'lbf':vec([1,1,-2,0,0]),
    'Pa':vec([1,1-2,-2,0,0]),
    'm':vec([0,1,0,0,0]),
    'ft':vec([0,1,0,0,0]),
    'm/s':vec([0,1,-1,0,0]),
    'g/m^3':vec([1,-3,0,0,0]),
    'kg/m^3':vec([1,-3,0,0,0]),
    'm/s^2':vec([0,1,-2,0,0]),
    'g_earth':vec([0,1,-2,0,0]),
    'kg':vec([1,0,0,0,0]),
    'lb_mass':vec([1,0,0,0,0]),
    'kg/s':vec([1,-1,0,0,0]),
    'W':vec([1,2,-3,0,0]),
    'hp':vec([1,2,-3,0,0]),
    'hz':vec([0,0,-1,0,0]),
    's':vec([0,0,1,0,0]),
    'm^3':vec([0,3,0,0,0]),
    'A':vec([0,0,0,1,0]),
    'C':vec([0,0,1,1,0]),
    'V':vec([1,2,-3,-1,0]),
    'J':vec([1,2,-2,0,0]),
    'K':vec([0,0,0,0,1]),
    'Ω':vec([1,2,-3,-2,0]),
    '':vec([0,0,0,0,0]),
    'H':vec([1,2,-2,-2,0]),
    'F':vec([-1,-2,4,2,0]),
    'heatflux':vec([0,-1,0,0,-1])+vec([1,2,-3,0,0]),#watts/meter*kelvin
    'conductivity':vec([0,-2,0,0,0])+vec([1,2,-3,0,0]),#watts/meter^2
    }

special_case_grams = ['g','g/m^3']

special_cases = {'g':.001,'g/m^3':.001,'lbf':4.44822,
    'hp':745.7,'ft':.3048,'lb_mass':0.453592,'g_earth':9.80665
}


def find_p(v):
    print('***************',v)
    best = 100000000000000
    prefix = ''
    # v = round(v,6)
    for p in list(prefixes):
        
        c = prefixes[p]
        vc = v/c
        
        #print('v',v,'c',c,'vc',vc,'best',best,'p',p)
        if vc<best and vc>=.99:
            prefix = p
            best = vc+0
            
    return(prefix)

def parallel(A,B):
    print (A,B)
    if ((A == 0) == (B==0)).all():
        print('true')
        print(A[A!=0])
        r = A[A!=0]/B[B!=0]
        if (r==r[0]).all():
            return True
    else:
        print('false')
    return False

def find_u(A):
    unit = ''
    for u in list(units):
        # print (A,units[u])
        if (A==units[u]).all():
            # print('true')
            unit = u
            break
    return (unit)

class o:
    def __init__(self,p,u):
        self.p = p
        self.u = u
        


def interms(a,b):
    A = np.transpose(np.stack((units[a[0]],units[a[1]],units[a[2]],units[a[3]])))
    X = inv(A).dot(b)
    
    for x in range(mlta):
        print(a[x]+'^'+'('+str(int(np.round(X[x])))+')',end='')
        if x<2:
            print(' * ',end='')
    
    return X
def interms_str(a,b):
    X = b
    # a = ['g','m','s','A','K']
    num = ''
    for x in range(mlta):
        if X[x] > 0:
            num += (a[x]+'*')*int(X[x])
    if num != '':
        if num[-1]=='*':
            num = num[0:-1]
    den = ''
    for x in range(mlta):
        if X[x] < 0:
            den += (a[x]+'*')*int(-X[x])
    if den != '':
        if den[-1]=='*':
            den = den[0:-1]
    if num == '':
        if den == '':
            return''
        else:
            return '1/'+den
    else:
        if den =='':
            return num
        else:
            return num+'/'+den

    return st

def multiply(pa,ua,pb,ub):
    prefix = find_p(prefixes[pa]*prefixes[pb])
    print(prefix)
    unit = interms(['kg','m','s'],units[ua]+units[ub])
    print( find_u(u))

def mult(L):
    p = 1
    u = vec([0,0,0])
    for l in L:
        p*=prefixes[l[0]]
        u+=units[l[1]]
    
    prefix = find_p(p)
    print()
    print(prefix,end='')
    # unit = interms(['kg','m','s'],u)
    print( find_u(u))


#multiply('k','kg','M','s')


def div(a,b):
    d=0
    print()
    v = 1
    if len(a)==3:
        #print ('initial',a[0]/b[0])
        d=1
        v = a[0]/b[0]
    
    
    p = prefixes[a[d+0]]
    u = units[a[d+1]]
    
    p/=prefixes[b[d+0]]
    u-=units[b[d+1]]
    
    show(v,p,u)

def show(v,p,u):
    
    unit = find_u(u)
    if find_u(u) == '':
        unit = interms_str(['g','m','s','A','K'],u)
        if len(unit)>0:
            if unit[0]=='g':
                gs = 0
                for un in unit:
                    if un == 'g':
                        gs+=1
                        
                p*= 1000**gs
                
            else:
                unit = interms_str(['kg','m','s','A','K'],u)
        ###########################################
    elif (u != units[find_u(u)]).any():
        unit = 'error'
    
    if unit in special_cases:
        p/=special_cases[unit]
    
    
    if unit == 'm^3':
        p = p**(1/3)
    
    # pre1 = prefixes[find_p(p)]
    pre1 = prefixes[find_p(p)]
    pre1 = p
    print('pre1',pre1)
    l10 = int(np.log10(abs(v)))
    print(l10,'l10')
    
    if unit == 'm^3':
        p*=10 ** (l10/3)
    else:
        p*= 10 ** l10
    
    print('v',v,'p',p)
    prefix = find_p(p)
    pre = prefixes[prefix]
    print('pre',pre)
    if unit == 'm^3':
        
        v *= (pre1/pre)**(3)
    else:
        v *= (pre1/pre)
    
    # print('v',v,'p',p)
    prefix = find_p(p)
    
    print(round(v,5),end = '  ')
    print(prefix,end='')
    # unit = interms(['kg','m','s'],u)
    print( unit)
    
    
    
    return[v,prefix,unit]


def m(L):
    v=1
    
    p = 1
    u = vec([0,0,0])
    for l in L:
        v*=l[0]
        p*=prefixes[l[1]]
        u+=units[l[2]]
    
    
    
    return show(v,p,u)

def m_d(L):
    v=1
    
    p = 1
    u = vec([0]*mlta)
    for l in L:
        
        if l[3] == 0:
            v*=l[0]
            p*=prefixes[l[1]]
            u+=units[l[2]]
            
            if l[2] in special_cases:
                p *= special_cases[l[2]]
            if l[2] == 'm^3':
                
                p*= prefixes[l[1]]*prefixes[l[1]]
            
        if l[3]==1:
            v/=l[0]
            p/=prefixes[l[1]]
            u-=units[l[2]]
            
            if l[2] in special_cases:
                p/=special_cases[l[2]]
            
            if l[2] == 'm^3':
                
                p/= prefixes[l[1]]*prefixes[l[1]]
    
    
    
    return show(v,p,u)

# m([
#     [1000,'k','m'],
#     [1000,'','m'],
#     [1,'','m'],
# ])

# m([[10,'','hz'],
# [10,'','hz'],
# [100,'c','m']])

from tkinter import *

from tkinter import ttk

# roo = Tk()
# frm1 = ttk.Frame(roo, padding=10)

# frm1.grid()



# ttk.Label(frm1, text="variable count").grid(column=0, row=0)

# variablecount = Text(frm1,height=1,width=6)
# variablecount.grid(column=1,row=0)
# variablecount.insert(INSERT, '1')
# def change_count():
#     # frm.destroy()
#     make_rows(int(variablecount.get(1.0, "end-1c")))

    
# ttk.Button(frm1, text="Apply", command=change_count).grid(column=2, row=0)

##########ROW 1
n=5
import sys
if len(sys.argv)>1:
    n = int(sys.argv[1])

root = Tk()
root.minsize(1520,300)
frm = ttk.Frame(root, padding=10)
frm.grid()
dpdwn_p=[None]*n
dpdwn_u=[None]*n

P = ['']*n

U = ['']*n
V = [None]*n
prin = [None]*n

printButton = [None]*n
dividecheck = [None]*n
divbool = [None]*n

# enter_equation = Text(frm,height=1,width=20)
# enter_equation.grid(column=1,row=0)
# enter_equation.insert(INSERT, '')
# enter_butt =Button(text="Add term", command=add_term)


for x in range(0,n):
    ttk.Label(frm, text="Term "+str(x)).grid(column=0, row=x+1)
    
    V[x] = Text(frm,height=1,width=6)
    V[x].grid(column=1,row=x+1)
    V[x].insert(INSERT, '1')
    
    P[x] = StringVar()
    dpdwn_p[x] = OptionMenu(frm,P[x],*list(prefixes))
    dpdwn_p[x].grid(column=2,row=1+x)
    U[x] = StringVar()
    U[x].set('')
    dpdwn_u[x] = OptionMenu(frm,U[x],*list(units))
    dpdwn_u[x].grid(column=3,row=1+x)
    
    divbool[x] = IntVar()
    dividecheck[x] = Checkbutton(frm,text = 'Divide?',variable = divbool[x],onvalue=1,offval=0)
    dividecheck[x].grid(column=4,row=1+x)

def get(x):
    inp = V[x].get(1.0, "end-1c")
    return[float(inp),P[x].get(),U[x].get(),divbool[x].get()]
def show_lab(x):
    
    print(get(x))
 
def show0():
    show_lab(0)
def show1():
    show_lab(1)
printButton[0] = ttk.Button(frm,
                        text = "Print",
                        command = show0).grid(column=5,row=1+0)
# printButton[1] = ttk.Button(frm,
#                         text = "Print",
#                         command = show1).grid(column=5,row=1+1)

answer_label = ttk.Label(frm,text='result')
answer_label.grid(column=3,row=1+n+1)

alt_answer = ttk.Label(frm,text='')
alt_answer.grid(column=4,row=1+n+1)

# form_denom = ttk.Label(frm,text='form denomenator')
# form_denom.grid(column=2,row=2+n+1)
formula_label = ttk.Label(frm,text='formula')
formula_label.grid(column=2,row=1+n+1)


def multipy_form():
    L = [get(x) for x in range(n)]
    # print(L)
    res = m_d(L)
    form = ''
    
    form_d = ''
    for l in L:
        print(l)
        
        if l[0] == 1.0 and l[1] == '' and l[2] == '' and l[3] == 0:
            pass
        else:
            if l[3] == 1:
                form_d+=str(round(l[0],5))
                form_d+=l[1]
                form_d+=l[2]
                form_d+= ' * '
                
                form+='('+str(round(l[0],5))+l[1]+l[2]+')^-1'
            
            else:
                form+=str(round(l[0],5))
                form+=l[1]
                form+=l[2]
            form+= ' * '
            
    form += ' = '
    
    
    
    
    
    # form_denom.config(text = form_d)
    formula_label.config(text = form)
    answer_label.config(text = str(round(res[0],5))+res[1]+res[2])
    
    # print('res',res[2])
    # if res[2] == 'm^3':
    
    #     alt_answer.config(text = 'also equals: ' + str(round(res[0]*1000,5))+res[1]+'Liter')
    # else:
    #     alt_answer.config(text = '')
    
    
def add_term():
    #import os
    import subprocess
    
    root.destroy()
    subprocess.call(['python3','prefixes.py',str(n+1)])
    
def remove_term():
    #import os
    import subprocess
    
    root.destroy()
    subprocess.call(['python3','prefixes.py',str(n-1)])


ttk.Button(frm, text="Add term", command=add_term    ).grid(column=4, row=1+1+n+1)
ttk.Button(frm, text="Remove term", command=remove_term    ).grid(column=3, row=1+1+n+1)


ttk.Button(frm, text="Multiply", command=multipy_form    ).grid(column=1, row=1+n+1)








# n+=1


ttk.Button(frm, text="Quit", command=root.destroy).grid(column=5, row=1+n+2+1)
root.mainloop()

            

# make_rows(3)


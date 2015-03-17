cd /Users/dmitry/gdisk/PHD/nla3d/
ls
cd custom/
ls
less model_lambda_20.txt
%less model_lambda_20.txt
np.genfromtxt?
np.genfromtxt("model_lambda_20.txt",skip_header = True)
M = np.genfromtxt("model_lambda_20.txt",skip_header = True)
shape(M)
his_values = []
his_weight = []

his_values = []
his_weight = []
for i in range(shape(M)[0]):
    lambda1 = M[i,2]
    lambda2 = M[i,3]
    lambda3 = M[i,4]
    J = lambda1*lambda2*lambda3
    lambda1 = J**(-1.0/3.0)*lambda1
    lambda2 = J**(-1.0/3.0)*lambda2
    lambda3 = J**(-1.0/3.0)*lambda3
    a = linspace(1.0,2.0,100)
    y1 = a**(1.0/2.0)-1.0/lamda2
    y2 = lambda1/lambda2 - a**(-3.0/2.0)
    plot(a,y1)
    plot(a,y2)
    raw_input("wait..")



i = 1000
    lambda1 = M[i,2]
    lambda2 = M[i,3]
    lambda3 = M[i,4]
    J = lambda1*lambda2*lambda3
    lambda1 = J**(-1.0/3.0)*lambda1
    lambda2 = J**(-1.0/3.0)*lambda2
    lambda3 = J**(-1.0/3.0)*lambda3
    print lambda1
    print lambda2
    print lambda3
    a = linspace(1.0,2.0,100)
    #y1 = a**(1.0/2.0)-1.0/lambda2
    #y2 = lambda1/lambda2 - a**(-3.0/2.0)
    y1 = a**(1.0/2.0)*lambda2-1
    y2 = a**(-3.0/2.0)*lambda1/lambda2-1
    plot(a,y1,"r-")
    a1 = max([lambda2**(-2), 1.0])
    a2 = max([(lambda2/lambda1)**(-2.0/3.0), 1.0])
    delta_a = abs(a1-a2)
    print "delta_a = " + str(delta_a)
    plot([a1],[0.0],'ro')
    plot([a2], [0.0], 'ro')
    plot(a,y2,"r--")
    grid()
    title("a,b,c range")
    xlabel("a,b,c")
    ylabel(">0")
    y1 = a**2*lambda2**(-2)-1
    y2 = a**(-3)*lambda1*lambda2**2-1
    plot(a,y1,"g-")
    plot(a,y2,"g--")
    b1 = max([lambda2, 1.0])
    b2 = max([(lambda1*lambda2**2)**(1.0/3.0),1.0])
    delta_b = abs(ab-ab)
    print "delta_b = " + str(delta_b)
    plot([b1],[0.0],'go')
    plot([b2],[0.0],'go')
    y1 = (a*lambda2/lambda1)**(-2.0/6.0)*lambda2-1
    c1 = 1.0
    c2 = lambda2**2*lambda1
    delta_c = abs(c1-c2)
    print "delta_c = " + str(delta_c)
    plot(a,y1,"b-")
    plot([c1],[0.0],'bo')
    plot([c2],[0.0],'bo')





al = []
bl = []
cl = []
wl = []
for i in range(shape(M)[0]):
    lambda1 = M[i,2]
    lambda2 = M[i,3]
    lambda3 = M[i,4]
    J = lambda1*lambda2*lambda3
    lambda1 = J**(-1.0/3.0)*lambda1
    lambda2 = J**(-1.0/3.0)*lambda2
    lambda3 = J**(-1.0/3.0)*lambda3
    I1 = lambda1**2+lambda2**2+lambda3**2
    wl.append(I1*M[i,1])
    #print lambda1
    #print lambda2
    #print lambda3
    a = linspace(1.0,2.0,100)
    #y1 = a**(1.0/2.0)-1.0/lambda2
    #y2 = lambda1/lambda2 - a**(-3.0/2.0)
    y1 = a**(1.0/2.0)*lambda2-1
    y2 = a**(-3.0/2.0)*lambda1/lambda2-1
    #plot(a,y1,"r-")
    a1 = max([lambda2**(-2), 1.0])
    a2 = max([(lambda2/lambda1)**(-2.0/3.0), 1.0])
    delta_a = abs(a1-a2)
    #print "delta_a = " + str(delta_a)
    #plot([a1],[0.0],'ro')
    #plot([a2], [0.0], 'ro')
    #plot(a,y2,"r--")
    #grid()
    #title("a,b,c range")
    #xlabel("a,b,c")
    #ylabel(">0")
    y1 = a**2*lambda2**(-2)-1
    y2 = a**(-3)*lambda1*lambda2**2-1
    #plot(a,y1,"g-")
    #plot(a,y2,"g--")
    b1 = max([lambda2, 1.0])
    b2 = max([(lambda1*lambda2**2)**(1.0/3.0),1.0])
    delta_b = abs(b1-b2)
    #print "delta_b = " + str(delta_b)
    #plot([b1],[0.0],'go')
    #plot([b2],[0.0],'go')
    y1 = (a*lambda2/lambda1)**(-2.0/6.0)*lambda2-1
    c1 = 1.0
    c2 = lambda2**2*lambda1
    delta_c = abs(c1-c2)
    #print "delta_c = " + str(delta_c)
    #plot(a,y1,"b-")
    #plot([c1],[0.0],'bo')
    #plot([c2],[0.0],'bo')
    al.append(delta_a)
    bl.append(delta_b)
    cl.append(delta_c)




ifile = open('newDump.txt','r')
#lowPt = 30
#ifile = open('dump_tt.txt','r')
lowPt = 40

zttInWindow = 0.
zttInWindowWt = 0.
zttInWindowGenPtMatch = 0.
zttInWindowGenPtMatchWt = 0.
zttAll = 0.
zttAllWt = 0.
zttAllGenPtMatch = 0.
zttAllGenPtMatchWt = 0.


for line in ifile :
    info = line.strip('\n').split(' ')

    # Make sure we have the right number of items (not a header row or something)
    if len(info) < 8 : continue
    # Only keep good ZTT gen matched events
    if info[7] == 'Gen_Match_Not5' : continue
    # Make sure the row doesn't have a non-numeric character for GenMass
    if info[0][0] not in ['0','1','2','3','4','5','6','7','8','9'] : continue


    if float(info[0]) > 60 and float(info[0]) < 120 :
        zttInWindow += 1
        zttInWindowWt += float(info[1])
        if (float(info[2]) > lowPt and float(info[3]) > 40) or (float(info[2]) > 40 and float(info[3]) > lowPt) :
            zttInWindowGenPtMatch += 1
            zttInWindowGenPtMatchWt += float(info[1])
            
    zttAll += 1
    zttAllWt += float(info[1])
    if (float(info[2]) > lowPt and float(info[3]) > 40) or (float(info[2]) > 40 and float(info[3]) > lowPt) :
        zttAllGenPtMatch += 1
        zttAllGenPtMatchWt += float(info[1])
        
print "zttInWindow: ",zttInWindow
print "zttInWindowWt: ",zttInWindowWt
print "zttInWindowGenPtMatch: ",zttInWindowGenPtMatch
print "zttInWindowGenPtMatchWt: ",zttInWindowGenPtMatchWt
print "zttAll: ",zttAll
print "zttAllWt: ",zttAllWt
print "zttAllGenPtMatch: ",zttAllGenPtMatch
print "zttAllGenPtMatchWt: ",zttAllGenPtMatchWt

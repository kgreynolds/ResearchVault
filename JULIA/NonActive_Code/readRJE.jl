
module readRJE #these functions are for reading regular input files

    function readDat(f)
        x = [];y=[]
        s = open(f) do file
            readlines(file)
        end
        for line in s
            spli = split(line)
            push!(x,parse(Float64,spli[1]));push!(y,parse(Float64,spli[2]))
        end
        return x,y
    end

    function determineDlm(lines)
        spacecount = 0; commacount = 0
        for l in lines
            if findfirst(" ",l) != nothing
                spacecount = 1 + spacecount
            elseif findfirst(",",l) != nothing
                commacount = 1 + commacount
            end
        end
        space = spacecount/length(lines)*100; comma = commacount/length(lines)*100
        println("dlm breakdown: $space% space delimited, $comma% comma delimited")
        if comma > space
            println("determined comma delimiting is more appropriate")
            return true
        else
            println("space delimiting is more appropriate")
            return false
        end
    end

    function readTxt(f) #this is designed to read any columned text file
        #specifically ones without decoration, but with headers/or not
        s = open(f) do file
            readlines(file)
        end

        #s = split.(s,",")#added this to make this play nicely with csvs
        if determineDlm(s)
            splt = split.(s,",")
        else
            splt = split.(s)
        end

        names = popfirst!(splt); noNames = false #assume headers
        out = popfirst!(splt) #get "seed"
        for s in splt
            out = hcat(out,s) #grow from "seed"
        end
        out = parse.(Float64,out) #parse the whole thing
        try #see if there are names
            names = parse.(Float64,names)
            out = hcat(names,out) #if there aren't, add those numbers back
            noNames = true
        catch #oops, there are headers made of words
            print("In fact, there are names: $names")
        end

        Out = [] #this is what I really plan on returning.
        for row in 1:size(out)[1] #iterates through the rows and collects them as Vecs
            push!(Out,out[row,:])
        end
        if noNames == true
            return Out
            print("did it")
        else
            return names,Out
        end
    end

    function readGr(file)
        lines = open(file) do file
            readlines(file)
        end
        startpoint = 0
        preSplt = split.(lines)
        #cycles throgh the lines that make the header of the file
        #can be modified in the future to pull out useful information
        for (i,splitline) in enumerate(preSplt)
            try
                if splitline == []
                    continue
                end
                attempt = parse.(Float64,splitline)
                startpoint = i
                #for (j,lin) in enumerate(lines[i:end])
                #end
                break
            catch
                #do literally nothing
                continue
            end
        end
        #printstyled("startpoint = $startpoint",color = :red)
        names = preSplt[startpoint-1][2:end]
        splt = preSplt[startpoint:end]
        out = popfirst!(splt) #get "seed"
        for s in splt
            out = hcat(out,s) #grow from "seed"
        end
        out = parse.(Float64,out) #parse the whole thing
        Out = [] #this is what I really plan on returning.
        for row in 1:size(out)[1] #iterates through the rows and collects the columns as Vecs
            push!(Out,out[row,:])
        end
        return names,Out #returns columns of data as Vecs with Headers
    end

    function readGrII(file)
        lines = open(file) do file
            readlines(file)
        end
        startpoint = 1
        #preSplt = split.(lines)
        if determineDlm(lines)
            preSplt = split.(lines,",")
        else
            preSplt = split.(lines)
        end
        #cycles throgh the lines that make the header of the file
        #can be modified in the future to pull out useful information
        for (i,splitline) in enumerate(preSplt)
            try
                if splitline == []
                    continue
                end
                attempt = parse.(Float64,splitline)
                startpoint = i
                print("found the startpoint!")
                try
                    for v in i:(i + trunc(Int,(length(lines))/20))
                        v0 = parse.(Float64,preSplt[v])
                        #printstyled("#$v worked: $v0",color = :green)
                    end
                catch
                    printstyled("lol, too soon \n ",color = :red)
                    continue
                end
                break #actually exit the loop having found what we were looking for
            catch
                print("not yet \n ")
                continue
            end

        end
        printstyled("startpoint = $startpoint",color = :red)
        metadata = lines[1:startpoint]
        #names = preSplt[startpoint - 1][2:end]
        splt = preSplt[startpoint:end]
        #printstyled(splt,color = :yellow)
        out = parse.(Float64,popfirst!(splt))
        printstyled(out,color = :red) #get "seed"
        for s in splt
            #printstyled(s,color = :blue)
            try
                out = vcat.(out,parse.(Float64,s)) #grow from "seed"
            catch
                print("$s: it seems there was a gap, or perhaps the end")
                continue
            end
        end
        #out = parse.(Float64,out)
        return metadata,out
        #return names,Out #returns columns of data as Vecs with Headers
    end

    function readCsv(file) #assumes column names and no other cells filled besides those with Float64 #s
        lines = open(file) do file
            readlines(file)
        end
        startpoint = 1
        preSplt = split.(lines,',')
        #cycles throgh the lines that make the header of the file
        #can be modified in the future to pull out useful information
        for (i,splitline) in enumerate(preSplt)
            try
                if splitline == []
                    continue
                end
                attempt = parse.(Float64,splitline)
                startpoint = i
                #for (j,lin) in enumerate(lines[i:end])
                #end
                break
            catch
                #do literally nothing
                continue
            end
        end
        #printstyled("startpoint = $startpoint",color = :red)
        names = preSplt[startpoint]
        splt = preSplt[startpoint+1:end]
        out = popfirst!(splt)[1:end-1] #get "seed", discarding empty terminal value
        for s in splt
            out = hcat(out,s[1:end-1]) #grow from "seed", again neglecting last value
        end
        out = parse.(Float64,out) #parse the whole thing
        Out = [] #this is what I really plan on returning.
        for row in 1:size(out)[1] #iterates through the rows and collects the columns as Vecs
            push!(Out,out[row,:])
        end
        return names,Out #returns columns of data as Vecs with Headers
    end

    function hashSQUID(file)
        startpoint = 0
        startpoint::Int
        datastartpoint = 0
        datastartpoint::Int

        lines = open(file) do file
            readlines(file)
        end

        for (i,line) in enumerate(lines) #this skips the bullshit and finds where
            #the data starts.
            if line == "[Data]"
                startpoint = i
            elseif line[1] == ','
                #println("found it")
                datastartpoint = i
                break
            else
                continue
            end
        end

        s = lines[startpoint+1:end]

        if determineDlm(s)
            splt = split.(s,",")
        else
            splt = split.(s)
        end

        names = popfirst!(splt); noNames = false #assume headers
        out = popfirst!(splt) #get "seed"
        for s in splt
            out = hcat(out,s) #grow from "seed"
        end

        try
            out = parse.(Float64,out) #parse the whole thing
        catch #empty fields don't parse well, gotta skip it
            println("contains missing values")
        end


        try #see if there are names
            names = parse.(Float64,names)
            out = hcat(names,out) #if there aren't, add those numbers back
            noNames = true
        catch #oops, there are headers made of words
            println("In fact, there are names")
        end

        Out = [] #this is what I really plan on returning.
        for row in 1:size(out)[1] #iterates through the rows and collects them as Vecs
            push!(Out,out[row,:])
        end
        if noNames == true
            return Out
            println("did it")
        else
            return names,Out
        end
        #return lines[startpoint:datastartpoint]
    end


    function siftParse(VecOfSubStrings)
        #oppotunistically parse. Some vectors will parse some won't.
        AnyVec = []
        for v in VecOfSubStrings
            try
                push!(AnyVec,parse.(Float64,v))
            catch
                push!(AnyVec,v)
            end
        end
        return AnyVec
    end

    function readSQUID(file, emptyBool)
        out = [[],[]]; goods = hashSQUID(file)
        goods = goods[1],siftParse(goods[2])
        #set things up for par
        if emptyBool
            for (i,g) in enumerate(goods[2])
                if typeof(g) == Vector{Float64}
                    push!(out[1],goods[1][i]);push!(out[2],g)
                    #keeps names and values of parsed data
                else
                    continue
                end
            end
        else
            out = goods # no tossing done
        end
        return out
    end

end  # module readRJE

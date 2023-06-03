#!/bin/bash
#
#
# 1. Sampler life times
#   - Start: AtlasBuilderActor Sent Sample message for sampling node
#   - SamplerActor received Sampling message for node
#   - SamplerActor sent sampling results for node
#   - Finish: SamplerActor received KILL message for node
# 2. Wait times in the queue
#   - Start: Adding to Queue, node
#   - Finish: AtlasBuilderActor Sent Sample message for sampling node
# 3. SamplerKillWaitTime: Time for each sampler between sending Done and receiving Kill
#   - SamplerActor about to send DONE message for node
#   - AtlasBuilderActor received DONE message from
#   - AtlasBuilderActor Sent Kill message to 
#   - SamplerActor received KILL message for node
#   - SamplerActor sent sampling results for node
# 4. RegionDiscoveryWaitTime: Time between new region sent to ABA and receiving boundary from ABA (match with CaleyPointID) 
#   - SamplerActor sent a new region to AB Actor from node
#   - Received New Region from
#   - SamplerActor received boundary message for
#   - AtlasBuilderActor sending boundary message to
# 5. Number of messages in the queue, every time there's a change in the queue (addition, deletion)
# 6. Number of samplers active, every time a new sampler is spawned or it returns. Grep "Number of active samplers:"
# 7. Wait time for witness writing
#

usage() {
    echo "sh extractLogs <Input_log_file> <Input_log_file> ..."
    echo "Creates out file: Outputfile_operation1 Outfile_operation2 ..."
}

outfile=""

for arg in "$@"
do
    case $arg in
        -o=*|--output=*)
        outfile="${arg#*=}"
        shift # Remove --output= from processing
        ;;
        -s=*|--select=*)
        select="${arg#*=}"
        shift # Remove --select= from processing
        ;;    
        #*)
        #OTHER_ARGUMENTS+=("$1")
        #shift # Remove generic argument from processing
        #;;
    esac
done

if [ $# -lt 1 ]; then
    echo "Missing input log file"
    usage
    exit 1
fi

# 1. SamplerLifeTime
SamplerLifeTime() {
    m1="AtlasBuilderActor Sent Sample message for sampling node"
    m2="SamplerActor received KILL message for node:"
    grep "$m1\|$m2" $1 | \
        sed "s|$m1|Sampler_R_Start|" | \
        sed "s|$m2|Sampler_R_Kill|" | \
        cut -d' ' -f 2,5,6 >> $outfile1
    echo "Appended ${outfile1}"
}

# 3. SamplerKillWaitTime
SamplerKillWaitTime() {
    m1="SamplerActor about to send DONE message for node:"
    #m2="AtlasBuilderActor received DONE message from:"
    #m3="AtlasBuilderActor Sent Kill message to node:"
    m4="SamplerActor received KILL message for node:"
    #m5="SamplerActor sent sampling results for node:"
    #grep "$m1\|$m2\|$m3\|$m4\|$m5" $1 | \
    grep "$m1\|$m4" $1 | \
        sed "s|$m1|Sampler_S_Done|" | \
        #sed "s|$m2|ABA_R_Done|" | \
        #sed "s|$m3|ABA_S_Kill|" | \
        sed "s|$m4|Sampler_R_Kill|" | \
        #sed "s|$m5|Sampler_S_Results|" | \
        cut -d' ' -f 2,5,6 >> $outfile3
    echo "Appended ${outfile3}"
}

# 4. RegionDiscoveryWaitTime
RegionDiscoveryWaitTime() {
    m1="SamplerActor sent a new region to AB Actor from node:"
    #m2="Received New Region from:"
    #m3="SamplerActor received boundary message for node:"
    m4="AtlasBuilderActor sending boundary message to:"
    #grep "$m1\|$m2\|$m3\|$m4" $1 | \
    grep "$m1\|$m4" $1 | \
        sed "s|$m1|Sampler_S_Region|" | \
        #sed "s|$m2|ABA_R_Region|" | \
        #sed "s|$m3|Sampler_R_Boundary|" | \
        sed "s|$m4|ABA_S_Boundary|" | \
        sed "s|CayleyPoint: ||" | \
        sed "s|Flip: ||" | \
        sed "s|BoundaryNode: ||" | \
        sed "s|ContactGraph: ||" | \
        sed "s|NumMessagesInABQ: [0-9]*||" | \
        cut -d' ' -f 2,5- >> $outfile4
    echo "Appended ${outfile4}"
}

for arg; do
    echo "Processing $arg"
    if [ $select -eq 0 ] || [ $select -eq 1 ]
    then 
        outfile1="${outfile}_part1"
        echo 'time_stamp message node' > $outfile1
        SamplerLifeTime $arg 
    fi
    if [ $select -eq 0 ] || [ $select -eq 2 ] 
    then 
        outfile3="${outfile}_part3"
        echo 'time_stamp message node' > $outfile3
        SamplerKillWaitTime $arg 
    fi
    if [ $select -eq 0 ] || [ $select -eq 3 ] 
    then 
        outfile4="${outfile}_part4"
        echo 'time_stamp message node graph cayley_point flip boundary' > $outfile4
        RegionDiscoveryWaitTime $arg 
    fi
done


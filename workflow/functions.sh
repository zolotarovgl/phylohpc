# Define the function to read and echo configuration variables
read_config() {
    local config_file="$1"

    if [[ ! -f "$config_file" ]]; then
        echo "Error: Config file '$config_file' not found!"
        return 1
    fi

    #echo "------------------"
    #echo -e "CONFIG_FILE: ${config_file}\n"

    # Read variables from the configuration file
    while IFS='=' read -r key value; do
        # Skip empty lines and comments
        if [[ -n "$key" && ! "$key" =~ ^# ]]; then
            eval "${key}='${value}'"
            # Echo the variable and its value
            #echo "$key=$value"
        fi
    done < "$config_file"

    #echo "------------------"
}

check_file_exists() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        echo "Warning: File '$file' does not exist."
        return 1
    else
        echo "$file exists"
        return 0
    fi
}


get_ram() {
    sacct --format=JobID,MaxRSS%9 -j "$1"
}


# check the status of the family
check_status_family() {
    local PREF FAMILY FILEPREF STATUS_S STATUS_C
    if [[ $# -eq 1 ]]; then
        IFS="." read -r PREF FAMILY <<< "$1"
    elif [[ $# -eq 2 ]]; then
        PREF=$1
        FAMILY=$2
    else
        echo -e "Usage: check_status PREF FAMILY\n   or: check_status PREF.FAMILY"
        return 1
    fi
    FILEPREF=${PREF}.${FAMILY}
    STATUS_S=$(check_file "$SEARCH_DIR/${FILEPREF}.genes.list")
    STATUS_C=$(check_file "$CLUSTER_DIR/${FILEPREF}_cluster.tsv")
    echo -e "#prefix\tsearch\tclustering"
    echo -e "$FILEPREF\t$STATUS_S\t$STATUS_C"
}

check_all_families() {
    if [[ $# -ne 1 ]]; then
        echo "Usage: check_all_families GENEFAM.csv"
        return 1
    fi
    local GENEFAM="$1"
    while read PREF FAMILY; do
        check_status_family "$PREF" "$FAMILY" | grep -v '#' | awk '{print $1"\t"$2$3}'
    done < <(awk '{print $NF"\t"$1}' "$GENEFAM")
}

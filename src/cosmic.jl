const data_dir   = "../data"
const result_dir = "../result"
const model_dir  = "../model"


function build_train()

end

function build_dir()
    for dir in (data_dir,result_dir,model_dir)
        if !isdir(dir)
            @info("buiding directory $data_dir")
            mkdir(dir)
        end
    end
end

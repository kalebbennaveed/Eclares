module GeneralUtils

# Rotation around Yaw
function Rot_Yaw(ϕ)
    Rz = [[cos(ϕ) ;; -sin(ϕ) ;; 0];
            [sin(ϕ) ;; cos(ϕ) ;; 0];
             [0   ;;  0  ;; 1]]
    return Rz
end

end

import { hash } from "bcrypt";
import { Response } from "express";

export const SHA256 = (msg: string) => {
    return hash(msg, 10);
};

const SECRET_KEY = "alles";

export const jwt = (header: any, payload: any) => {
    return SHA256(String(header) + "." + String(payload) + "." + SECRET_KEY);
};

export const throwError = (stack: string, status: number, msg: string)=>{
    throw {
        stack: stack,
        status: status,
        message: msg
    }
}

export const errorResponce = (err: any, res: Response): void => {
    console.error(err.stack);
    res.status(err.status || 500).json({
        error: {
            message: err.message || "Internal Server Error",
        },
    });
};

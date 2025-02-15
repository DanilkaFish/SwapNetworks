import express from "express";
// const express = require('express');
import { register, login, refreshToken, logout } from "../controller/user.controller";
import { auth } from "../middleware/auth.middleware";

export const authRouter = express.Router();

authRouter.post("/register", register);
authRouter.post("/login", login);
authRouter.post("/refresh-token", refreshToken);
authRouter.post("logout", auth, logout)
